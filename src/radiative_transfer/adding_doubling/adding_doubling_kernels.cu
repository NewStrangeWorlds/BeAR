// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include "adding_doubling.h"

#include "../../CUDA_kernels/data_management_kernels.h"
#include "../../additional/exceptions.h"

#include <cuda_solver.cuh>

#include <thrust/device_ptr.h>
#include <thrust/extrema.h>

#include <string>


namespace bear{


namespace {

float maxOnDevice(const float* data, const size_t nb_elements)
{
  thrust::device_ptr<const float> ptr(data);
  return *thrust::max_element(ptr, ptr + nb_elements);
}

}


AddingDoubling::~AddingDoubling()
{
  if (phase_moments_dev != nullptr) deleteFromDevice(phase_moments_dev);

  if (gpu_workspace_ptr != nullptr)
    delete static_cast<adrt::cuda::SolverWorkspaceGPU*>(gpu_workspace_ptr);
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

  // Gas phase function: Rayleigh moments (chi_0 = 1, chi_2 = 0.1), shared across
  // wavelengths. Cloud scattering (Henyey-Greenstein) is mixed in by the solver.
  if (phase_moments_dev == nullptr)
  {
    std::vector<float> phase_moments_host(nb_layers * nb_phase_moments, 0.0f);
    for (size_t l = 0; l < nb_layers; ++l)
    {
      phase_moments_host[l * nb_phase_moments] = 1.0f;
      if (nb_phase_moments > 2) phase_moments_host[l * nb_phase_moments + 2] = 0.1f;
    }

    allocateOnDevice(phase_moments_dev, nb_layers * nb_phase_moments);
    moveToDevice(phase_moments_dev, phase_moments_host);
  }

  // The batched solver for more than 16 quadrature points only supports phase
  // moments shared across wavelengths and can therefore not treat scattering clouds
  if (nb_quadrature > 16 && cloud_single_scattering_dev != nullptr
      && maxOnDevice(cloud_single_scattering_dev, nb_layers * nb_spectral_points) > 0.0f)
  {
    std::string error_message =
      "Adding-doubling RT on the GPU supports scattering clouds only for up to 16 quadrature points\n";
    throw InvalidInput(std::string ("AddingDoubling::calcSpectrumGPU"), error_message);
  }

  // Configure the batched solver
  adrt::cuda::BatchConfig bcfg;
  bcfg.num_wavenumbers = nb_spectral_points;
  bcfg.num_layers = nb_layers;
  bcfg.num_quadrature = nb_quadrature;
  bcfg.num_moments_max = nb_phase_moments;
  bcfg.use_delta_m = use_delta_m;
  bcfg.surface_albedo = 0.0;
  bcfg.solar_flux = 0.0;
  bcfg.solar_mu = 0.5;
  bcfg.use_thermal_emission = true;
  bcfg.top_temperature = 0.0;  // cold space: no thermal downwelling at TOA
  bcfg.spectrum_scaling = spectrum_scaling;

  // Pass BeAR's raw data directly to the solver
  adrt::cuda::RawDeviceData data;
  data.absorption_coeff = absorption_coeff_dev;
  data.scattering_coeff = scattering_coeff_dev;
  data.altitude = atmosphere.altitude_dev;
  data.temperature = atmosphere.temperature_dev;
  data.wavenumber = spectral_grid->wavenumber_list_gpu;
  data.phase_moments = phase_moments_dev;
  data.phase_moments_shared = true;
  // clouds enter as the particle component; the library mixes them with the gas
  data.particle_optical_depth = cloud_optical_depth_dev;
  data.particle_single_scattering = cloud_single_scattering_dev;
  data.particle_asym_param = cloud_asym_param_dev;
  data.flux_up = model_spectrum_dev;

  if (nb_quadrature <= 16)
  {
    adrt::cuda::solveBatchFromCoefficients(bcfg, data);
  }
  else
  {
    if (gpu_workspace_ptr == nullptr)
      gpu_workspace_ptr = new adrt::cuda::SolverWorkspaceGPU();

    auto& ws = *static_cast<adrt::cuda::SolverWorkspaceGPU*>(gpu_workspace_ptr);
    ws.allocate(nb_spectral_points, nb_layers);
    adrt::cuda::solveBatchFromCoefficients(bcfg, data, ws);
  }
}


}
