// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <vector>
#include <cmath>

#include "discrete_ordinate_adding.h"
#include "dda_layer.h"

#include "../../forward_model/atmosphere/atmosphere.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../CUDA_kernels/error_check.h"
#include "../../CUDA_kernels/planck_function.h"


namespace bear{


//one thread per spectral point; the temperature and altitude profiles are staged in shared memory
//the layers are added from the bottom up, the layer algebra is shared with the CPU version (dda_layer.h)
//and evaluated in single precision: its cancellation-free form keeps the result within a few 1e-6
//of the double-precision CPU solution, also for thin layers and nearly conservative scattering
__global__
void discreteOrdinateAddingDev(
  float* __restrict__ model_spectrum_gpu,
  const float* __restrict__ absorption_coeff_dev,
  const float* __restrict__ scattering_coeff_dev,
  const float* __restrict__ cloud_optical_depth_dev,
  const float* __restrict__ cloud_single_scattering_dev,
  const float* __restrict__ cloud_asym_param_dev,
  const double* __restrict__ wavenumber_list_dev,
  const float* __restrict__ temperature_dev,
  const float* __restrict__ altitude_dev,
  const bool use_delta_m,
  const double spectrum_scaling,
  const int nb_spectral_points,
  const int nb_grid_points)
{
  extern __shared__ float shared_data[];
  float* s_altitude = shared_data;
  float* s_temperature = &shared_data[nb_grid_points];

  for (int i = threadIdx.x; i < nb_grid_points; i += blockDim.x)
  {
    s_altitude[i] = altitude_dev[i];
    s_temperature[i] = temperature_dev[i];
  }

  __syncthreads();

  const int nb_layers = nb_grid_points - 1;


  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_spectral_points; tid += blockDim.x * gridDim.x)
  {
    const float wavenumber = static_cast<float>(wavenumber_list_dev[tid]);
    const float wavenumber_cube = wavenumber * wavenumber * wavenumber;

    //lower boundary: black surface emitting at the bottom temperature
    float planck_bot = planckFunctionAccurate(s_temperature[0], wavenumber_cube, wavenumber);

    dda::Mat2<float> reflection = {0.0f, 0.0f, 0.0f, 0.0f};
    dda::Vec2<float> intensity = {planck_bot, planck_bot};

    for (int i = 0; i < nb_layers; ++i)
    {
      const float planck_top = planckFunctionAccurate(s_temperature[i+1], wavenumber_cube, wavenumber);

      //optical properties of layer i from the level-based gas coefficients and the layer-based cloud properties
      const int idx0 = i * nb_spectral_points + tid;
      const int idx1 = (i + 1) * nb_spectral_points + tid;

      const float dz = s_altitude[i+1] - s_altitude[i];

      const float absorption_depth = dz * (absorption_coeff_dev[idx1] + absorption_coeff_dev[idx0]) * 0.5f;
      const float gas_scattering_depth = dz * (scattering_coeff_dev[idx1] + scattering_coeff_dev[idx0]) * 0.5f;

      float tau = absorption_depth + gas_scattering_depth;
      float cloud_scattering_depth = 0.0f;
      float cloud_absorption_depth = 0.0f;
      float g_cloud = 0.0f;

      if (cloud_optical_depth_dev != nullptr)
      {
        const float cloud_tau = cloud_optical_depth_dev[idx0];
        tau += cloud_tau;
        cloud_absorption_depth = cloud_tau;

        if (cloud_single_scattering_dev != nullptr)
        {
          const float cloud_ssa = cloud_single_scattering_dev[idx0];
          cloud_scattering_depth = cloud_tau * cloud_ssa;
          cloud_absorption_depth = cloud_tau * (1.0f - cloud_ssa);
        }

        if (cloud_asym_param_dev != nullptr)
          g_cloud = cloud_asym_param_dev[idx0];
      }

      if (tau < 0.0f) tau = 0.0f;

      const float scattering_depth = gas_scattering_depth + cloud_scattering_depth;

      float omega = 0.0f;
      float one_minus_omega = 1.0f;

      if (tau > 0.0f)
      {
        omega = fminf(scattering_depth / tau, 1.0f);
        one_minus_omega = fminf(fmaxf((absorption_depth + cloud_absorption_depth) / tau, 0.0f), 1.0f);
      }

      const float cloud_fraction = (scattering_depth > 0.0f) ? cloud_scattering_depth / scattering_depth : 0.0f;

      float chi[dda::nb_moments];

      dda::scaledOpticalProperties(tau, omega, one_minus_omega, cloud_fraction, g_cloud, use_delta_m, chi);

      const dda::Layer<float> layer = dda::layerSolution(tau, omega, one_minus_omega, chi, planck_top, planck_bot);

      dda::addLayer(reflection, intensity, layer);

      planck_bot = planck_top;
    }

    model_spectrum_gpu[tid] = static_cast<float>(dda::upwardFlux(intensity) * 1e-3 * spectrum_scaling);
  }
}



void DiscreteOrdinateAdding::calcSpectrumGPU(
  const Atmosphere& atmosphere,
  float* absorption_coeff_dev,
  float* scattering_coeff_dev,
  float* cloud_optical_depth_dev,
  float* cloud_single_scattering_dev,
  float* cloud_asym_param_dev,
  const double spectrum_scaling,
  float* model_spectrum_dev)
{
  const size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  discreteOrdinateAddingDev<<<blocks, threads, 2 * nb_grid_points * sizeof(float)>>>(
    model_spectrum_dev,
    absorption_coeff_dev,
    scattering_coeff_dev,
    cloud_optical_depth_dev,
    cloud_single_scattering_dev,
    cloud_asym_param_dev,
    spectral_grid->wavenumber_list_gpu,
    atmosphere.temperature_dev,
    atmosphere.altitude_dev,
    use_delta_m,
    spectrum_scaling,
    static_cast<int>(nb_spectral_points),
    static_cast<int>(nb_grid_points));

  CUDA_CHECK_AFTER_KERNEL();
}


}
