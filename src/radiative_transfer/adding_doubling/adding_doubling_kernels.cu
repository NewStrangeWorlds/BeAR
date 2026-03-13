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

#include <cuda_solver.cuh>


namespace bear{


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

  // Allocate isotropic phase moments on first call
  if (phase_moments_dev == nullptr)
  {
    nb_phase_moments = 2 * nb_quadrature;

    std::vector<float> phase_moments_host(nb_layers * nb_phase_moments, 0.0f);
    for (size_t l = 0; l < nb_layers; ++l)
      phase_moments_host[l * nb_phase_moments] = 1.0f;

    allocateOnDevice(phase_moments_dev, nb_layers * nb_phase_moments);
    moveToDevice(phase_moments_dev, phase_moments_host);
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
  bcfg.use_thermal_emission = true;
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
  data.cloud_optical_depth = cloud_optical_depth_dev;
  data.flux_up = model_spectrum_dev;

  if (nb_quadrature <= 8)
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
