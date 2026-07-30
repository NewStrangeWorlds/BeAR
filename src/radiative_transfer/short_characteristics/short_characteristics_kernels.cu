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


#include <iostream>
#include <vector>
#include "math.h"
#include <stdio.h>
#include <new>

#include "short_characteristics.h"

#include "../../forward_model/atmosphere/atmosphere.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../additional/physical_const.h"
#include "../../CUDA_kernels/error_check.h"
#include "../../CUDA_kernels/planck_function.h"


namespace bear{

//solves the radiative transfer equation with the short characteristic method
//uses two angles, distributed according to a Gaussian quadrature scheme
//the temperature and altitude profiles are staged in shared memory
__global__
void shortCharacteristicsDev_Shared(
  float* __restrict__ model_spectrum_gpu,
  const float* __restrict__ absorption_coeff_dev,
  const double* __restrict__ wavenumber_list_dev,
  const float* __restrict__ cloud_optical_depth_dev,
  const float* __restrict__ temperature_dev,
  const float* __restrict__ vertical_grid_dev,
  const double spectrum_scaling,
  const int nb_spectral_points,
  const int nb_grid_points)
{
  extern __shared__ float shared_data[];
  float* s_vertical_grid = shared_data;
  float* s_temperature = &shared_data[nb_grid_points];

  // Collaborative loading: all threads help load the atmospheric profile
  for (int i = threadIdx.x; i < nb_grid_points; i += blockDim.x)
  {
    s_vertical_grid[i] = vertical_grid_dev[i];
    s_temperature[i] = temperature_dev[i];
  }

  __syncthreads(); // Ensure the profile is fully loaded before proceeding

  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_spectral_points; tid += blockDim.x * gridDim.x)
  {
    const float mu1 = 0.211324865405187;
    const float mu2 = 0.788675134594813;
    const float wavenumber = wavenumber_list_dev[tid];
    const float wn_cube = wavenumber * wavenumber * wavenumber;

    // Use shared memory for initial boundary condition
    float p_curr = planckFunction(s_temperature[0], wn_cube, wavenumber);
    double intensity_mu1 = p_curr;
    double intensity_mu2 = p_curr;

    float z_curr = s_vertical_grid[0];
    float abs_curr = absorption_coeff_dev[tid];

    for (int i = 0; i < nb_grid_points - 1; ++i)
    {
      const float z_next = s_vertical_grid[i+1];
      const float p_next = planckFunction(s_temperature[i+1], wn_cube, wavenumber);

      // Global memory load (absorption is unique per spectral point/thread);
      // the current layer's value is carried over from the previous iteration
      const float abs_next = absorption_coeff_dev[(i + 1) * nb_spectral_points + tid];

      float tau_layer = (z_next - z_curr) * (abs_next + abs_curr) * 0.5;

      if (cloud_optical_depth_dev != nullptr)
        tau_layer += cloud_optical_depth_dev[i * nb_spectral_points + tid];

      if (tau_layer > 1e-12f)
      {
        // Mu 1 path
        float d1 = tau_layer / mu1;
        float att1 = __expf(-d1);
        float term1 = expm1f(-d1) / d1;
        intensity_mu1 = intensity_mu1 * att1 + (1.0f + term1) * p_next + (-att1 - term1) * p_curr;

        // Mu 2 path
        float d2 = tau_layer / mu2;
        float att2 = __expf(-d2);
        float term2 = expm1f(-d2) / d2;
        intensity_mu2 = intensity_mu2 * att2 + (1.0f + term2) * p_next + (-att2 - term2) * p_curr;
      }

      z_curr = z_next;
      p_curr = p_next;
      abs_curr = abs_next;
    }

    const double final_const = 3.141592653589793 * 1e-3 * spectrum_scaling;
    model_spectrum_gpu[tid] = static_cast<float>(final_const * (intensity_mu1 * mu1 + intensity_mu2 * mu2));
  }
}



void ShortCharacteristics::calcSpectrumGPU(
  const Atmosphere& atmosphere,
  float* absorption_coeff_dev,
  float* scattering_coeff_dev,
  float* cloud_optical_depth_dev,
  float* cloud_single_scattering_dev,
  float* cloud_asym_param_dev,
  const double spectrum_scaling,
  float* model_spectrum_dev)
{
  size_t nb_grid_points = atmosphere.temperature.size();
  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;


  shortCharacteristicsDev_Shared<<<blocks,threads, 2 * nb_grid_points * sizeof(float)>>>(
    model_spectrum_dev,
    absorption_coeff_dev,
    spectral_grid->wavenumber_list_gpu,
    cloud_optical_depth_dev,
    atmosphere.temperature_dev,
    atmosphere.altitude_dev,
    spectrum_scaling,
    nb_spectral_points,
    nb_grid_points);

  CUDA_CHECK_AFTER_KERNEL();
}


}
