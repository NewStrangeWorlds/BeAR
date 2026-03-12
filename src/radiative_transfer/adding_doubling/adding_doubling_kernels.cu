/*
* This file is part of the BeAR code (https://github.com/newstrangeworlds/BeAR).
* Copyright (C) 2024 Daniel Kitzmann
*
* BeAR is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* BeAR is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* BeAR directory under <LICENSE>. If not, see
* <http://www.gnu.org/licenses/>.
*/


#include "adding_doubling.h"

#include "../../CUDA_kernels/data_management_kernels.h"
#include "../../CUDA_kernels/planck_function.h"
#include "../../CUDA_kernels/error_check.h"

#include <cuda_solver.cuh>


namespace bear{


// Compute optical depth and single scattering albedo per layer from
// BeAR's absorption and scattering coefficients at grid points.
//
// BeAR layout:  coeff[level * nwav + wav]  (level-major, level 0 = BOA)
// ADRT layout:  delta_tau[wav * nlay + layer]  (wavenumber-major, layer 0 = TOA)
//
// The ADRT CUDA solver has no index_from_bottom option, so we reverse
// the layer ordering here: ADRT layer 0 (top) = BeAR layer nlay-1.
//
// optical_depth = dz * (ext[j+1] + ext[j]) / 2 + cloud_tau
// where ext = absorption + scattering (extinction coefficient)
// SSA = scattering_depth / optical_depth
__global__ void computeOpticalPropertiesKernel(
  float* __restrict__ delta_tau_out,
  float* __restrict__ ssa_out,
  const float* __restrict__ absorption_coeff_dev,
  const float* __restrict__ scattering_coeff_dev,
  const float* __restrict__ cloud_optical_depth_dev,
  const float* __restrict__ altitude_dev,
  const int nb_spectral_points,
  const int nb_grid_points)
{
  const int nb_layers = nb_grid_points - 1;

  for (int tid = blockIdx.x * blockDim.x + threadIdx.x;
       tid < nb_spectral_points * nb_layers;
       tid += blockDim.x * gridDim.x)
  {
    const int wav = tid / nb_layers;
    const int bear_layer = tid % nb_layers;

    const float dz = altitude_dev[bear_layer + 1] - altitude_dev[bear_layer];

    // Extinction = absorption + scattering at each grid point
    const float ext_bot = absorption_coeff_dev[bear_layer * nb_spectral_points + wav]
                        + scattering_coeff_dev[bear_layer * nb_spectral_points + wav];
    const float ext_top = absorption_coeff_dev[(bear_layer + 1) * nb_spectral_points + wav]
                        + scattering_coeff_dev[(bear_layer + 1) * nb_spectral_points + wav];

    // Scattering contribution to layer optical depth
    const float scat_bot = scattering_coeff_dev[bear_layer * nb_spectral_points + wav];
    const float scat_top = scattering_coeff_dev[(bear_layer + 1) * nb_spectral_points + wav];
    const float scat_depth = dz * (scat_top + scat_bot) * 0.5f;

    float optical_depth = dz * (ext_top + ext_bot) * 0.5f;

    if (cloud_optical_depth_dev != nullptr)
      optical_depth += cloud_optical_depth_dev[bear_layer * nb_spectral_points + wav];

    if (optical_depth < 0.0f) optical_depth = 0.0f;

    // Reverse layer order: ADRT layer 0 = top = BeAR layer nlay-1
    const int adrt_layer = nb_layers - 1 - bear_layer;
    const int out_idx = wav * nb_layers + adrt_layer;
    delta_tau_out[out_idx] = optical_depth;
    ssa_out[out_idx] = (optical_depth > 0.0f) ? scat_depth / optical_depth : 0.0f;
  }
}


// Compute Planck levels for each wavenumber and atmospheric level.
//
// BeAR:  temperature_dev[level]  (level 0 = BOA)
// ADRT:  planck_levels[wav * nlev + level]  (level 0 = TOA)
//
// Level ordering is reversed to match the ADRT convention.
__global__ void computePlanckLevelsKernel(
  float* __restrict__ planck_levels_out,
  const float* __restrict__ temperature_dev,
  const double* __restrict__ wavenumber_list_dev,
  const int nb_spectral_points,
  const int nb_grid_points)
{
  for (int wav = blockIdx.x * blockDim.x + threadIdx.x;
       wav < nb_spectral_points;
       wav += blockDim.x * gridDim.x)
  {
    const float wn = static_cast<float>(wavenumber_list_dev[wav]);
    const float wn_cube = wn * wn * wn;

    for (int adrt_level = 0; adrt_level < nb_grid_points; ++adrt_level)
    {
      // Reverse level order: ADRT level 0 = TOA = BeAR level nlev-1
      const int bear_level = nb_grid_points - 1 - adrt_level;

      planck_levels_out[wav * nb_grid_points + adrt_level] = planckFunction(
        temperature_dev[bear_level], wn_cube, wn);
    }
  }
}


// Apply spectrum scaling factor to the output flux
__global__ void applyScalingKernel(
  float* __restrict__ spectrum,
  const float scaling,
  const int nb_points)
{
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x;
       tid < nb_points;
       tid += blockDim.x * gridDim.x)
  {
    spectrum[tid] *= scaling;
  }
}


void AddingDoubling::calcSpectrumGPU(
  const Atmosphere& atmosphere,
  float* absorption_coeff_dev,
  float* scattering_coeff_dev,
  float* cloud_optical_depth_dev,
  float* cloud_single_scattering_dev,
  float* cloud_asym_param_dev,
  const double spectrum_scaling,
  float* model_spectrum_dev)
{
  const size_t nb_layers = nb_grid_points - 1;
  const size_t nb_spectral_points = spectral_grid->nbSpectralPoints();
  const size_t nwav_nlay = nb_spectral_points * nb_layers;
  const size_t nwav_nlev = nb_spectral_points * nb_grid_points;

  // Allocate persistent device buffers on first call
  if (gpu_buffer_size == 0)
  {
    gpu_buffer_size = nb_spectral_points;

    // Phase moments: isotropic, shared across all wavenumbers
    nb_phase_moments = 2 * nb_quadrature;

    std::vector<float> phase_moments_host(nb_layers * nb_phase_moments, 0.0f);
    for (size_t l = 0; l < nb_layers; ++l)
      phase_moments_host[l * nb_phase_moments] = 1.0f;

    allocateOnDevice(phase_moments_dev, nb_layers * nb_phase_moments);
    moveToDevice(phase_moments_dev, phase_moments_host);

    // Working buffers for converted optical properties and Planck levels
    allocateOnDevice(delta_tau_dev, nwav_nlay);
    allocateOnDevice(ssa_dev, nwav_nlay);
    allocateOnDevice(planck_levels_dev, nwav_nlev);
  }

  const int threads = 256;

  // Compute optical depth and SSA from absorption/scattering coefficients
  {
    int blocks = (nwav_nlay + threads - 1) / threads;
    computeOpticalPropertiesKernel<<<blocks, threads>>>(
      delta_tau_dev,
      ssa_dev,
      absorption_coeff_dev,
      scattering_coeff_dev,
      cloud_optical_depth_dev,
      atmosphere.altitude_dev,
      nb_spectral_points,
      nb_grid_points);
    CUDA_CHECK_AFTER_KERNEL();
  }

  // Compute Planck function at each level for each wavenumber
  {
    int blocks = (nb_spectral_points + threads - 1) / threads;
    computePlanckLevelsKernel<<<blocks, threads>>>(
      planck_levels_dev,
      atmosphere.temperature_dev,
      spectral_grid->wavenumber_list_gpu,
      nb_spectral_points,
      nb_grid_points);
    CUDA_CHECK_AFTER_KERNEL();
  }

  // Configure the batched solver
  adrt::cuda::BatchConfig bcfg;
  bcfg.num_wavenumbers = nb_spectral_points;
  bcfg.num_layers = nb_layers;
  bcfg.num_quadrature = nb_quadrature;
  bcfg.num_moments_max = nb_phase_moments;
  bcfg.surface_albedo = 0.0;
  bcfg.solar_flux = 0.0;
  bcfg.solar_mu = 0.5;
  bcfg.use_thermal_emission = false;

  adrt::cuda::DeviceData data;
  data.delta_tau = delta_tau_dev;
  data.single_scat_albedo = ssa_dev;
  data.phase_moments = phase_moments_dev;
  data.phase_moments_shared = true;
  data.planck_levels = planck_levels_dev;
  data.flux_up = model_spectrum_dev;
  data.flux_down = nullptr;
  data.flux_direct = nullptr;

  adrt::cuda::solveBatch(bcfg, data);

  // Convert cgs (erg s-1 cm-2) to SI (W m-2) and apply spectrum scaling
  {
    const float scaling = static_cast<float>(1e-3 * spectrum_scaling);
    int blocks = (nb_spectral_points + threads - 1) / threads;
    applyScalingKernel<<<blocks, threads>>>(
      model_spectrum_dev,
      scaling,
      nb_spectral_points);
    CUDA_CHECK_AFTER_KERNEL();
  }
}


}
