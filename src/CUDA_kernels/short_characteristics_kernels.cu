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

#include "../radiative_transfer/short_characteristics.h"

#include "../forward_model/atmosphere/atmosphere.h"
#include "../spectral_grid/spectral_grid.h"
#include "../additional/physical_const.h"
#include "error_check.h"
#include "planck_function.h"


namespace bear{

//solves the radiative transfer equation with the short characteristic method
//uses two angles, distributed according to a Gaussian quadrature scheme
__global__ void shortCharacteristicsDevOld(
  double* model_spectrum_gpu,
  const float* absorption_coeff_dev,
  const double* wavenumber_list_dev,
  const float* cloud_optical_depth_dev,
  const float* temperature_dev,
  const float* vertical_grid_dev,
  const double spectrum_scaling,
  const int nb_spectral_points,
  const int nb_grid_points)
{
  //tid is the wavenumber index
  //int tid = blockIdx.x * blockDim.x + threadIdx.x;
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_spectral_points; tid += blockDim.x * gridDim.x)
  {
    constexpr double gauss_nodes[2] = {0.211324865405187, 0.788675134594813};
    constexpr double gauss_weights[2] = {0.5, 0.5};

    //inner boundary condition
    double intensity_mu1 = planckFunction(temperature_dev[0], wavenumber_list_dev[tid]);
    double intensity_mu2 = intensity_mu1;

    for (int i=0; i<nb_grid_points-1; ++i)
    { 
      const double planck_function1 = planckFunction(temperature_dev[i], wavenumber_list_dev[tid]);
      const double planck_function2 = planckFunction(temperature_dev[i+1], wavenumber_list_dev[tid]);

      //optical depth includes the molecular absorption coefficients as well as the one of the cloud
      double optical_depth_layer = 
        (static_cast<double>(vertical_grid_dev[i+1]) - static_cast<double>(vertical_grid_dev[i])) 
        * ( static_cast<double>(absorption_coeff_dev[(i+1)*nb_spectral_points + tid]) + static_cast<double>(absorption_coeff_dev[i*nb_spectral_points + tid]))/2.;
      
        if (cloud_optical_depth_dev != nullptr)
        optical_depth_layer += cloud_optical_depth_dev[i*nb_spectral_points + tid];

      if (optical_depth_layer == 0)
        continue;

      //Gauss angle 1
      const double delta1 = optical_depth_layer/gauss_nodes[0];
      const double attenuation_factor1 = exp(-delta1);

      intensity_mu1 = intensity_mu1 * attenuation_factor1;

      const double beta1 = 1.0 + (attenuation_factor1 - 1.0)/delta1;
      const double gamma1 = -attenuation_factor1 - (attenuation_factor1 - 1.0)/delta1;

      intensity_mu1 += beta1 * planck_function2 + gamma1 * planck_function1;


      //Gauss angle 2
      const double delta2 = optical_depth_layer/gauss_nodes[1];
      const double attenuation_factor2 = exp(-delta2);

      intensity_mu2 = intensity_mu2 * attenuation_factor2;

      const double beta2 = 1.0 + (attenuation_factor2 - 1.0)/delta2;
      const double gamma2 = -attenuation_factor2 - (attenuation_factor2 - 1.0)/delta2;

      intensity_mu2 += beta2 * planck_function2 + gamma2 * planck_function1;
    }

    //and integration with the corresponding Gauss-Legendre quadrature weights
    model_spectrum_gpu[tid] = 
    2.0 * constants::pi 
    * (intensity_mu1 * gauss_nodes[0] * gauss_weights[0] 
     + intensity_mu2 * gauss_nodes[1] * gauss_weights[1]) * 1e-3 * spectrum_scaling; //in W m-2 cm-1
  }
}


__global__ 
void shortCharacteristicsDev(
  double* __restrict__ model_spectrum_gpu,
  const float* __restrict__ absorption_coeff_dev,
  const double* __restrict__ wavenumber_list_dev,
  const float* __restrict__ cloud_optical_depth_dev,
  const float* __restrict__ temperature_dev,
  const float* __restrict__ vertical_grid_dev,
  const double spectrum_scaling,
  const int nb_spectral_points,
  const int nb_grid_points)
{
  // Grid-stride loop for flexibility with block sizes
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_spectral_points; tid += blockDim.x * gridDim.x)
  {
    // Use constants for Gauss nodes/weights
    const double mu1 = 0.211324865405187;
    const double mu2 = 0.788675134594813;
    
    const double wavenumber = wavenumber_list_dev[tid];
    
    // Initial Boundary Condition
    // Pre-calculate and cache Planck values to avoid redundant calls
    double p_curr = planckFunction(temperature_dev[0], wavenumber);
    double intensity_mu1 = p_curr;
    double intensity_mu2 = p_curr;

    // Cache altitude for the current layer
    double z_curr = vertical_grid_dev[0];

    for (int i = 0; i < nb_grid_points - 1; ++i)
    {
      // Fetch next layer data
      const float z_next = vertical_grid_dev[i+1];
      const float p_next = planckFunction(temperature_dev[i+1], wavenumber);
      
      // Optical depth calculation
      // Combining memory loads: fetch absorption for next layer
      float abs_curr = absorption_coeff_dev[i * nb_spectral_points + tid];
      float abs_next = absorption_coeff_dev[(i + 1) * nb_spectral_points + tid];
      
      float tau_layer = (z_next - z_curr) * (abs_next + abs_curr) * 0.5;

      if (cloud_optical_depth_dev != nullptr)
        tau_layer += cloud_optical_depth_dev[i * nb_spectral_points + tid];

      // Early exit if layer is transparent
      if (tau_layer > 1e-12) 
      {
        // Calculation for Mu 1
        double d1 = tau_layer / mu1;
        double att1 = exp(-d1);
        double inv_d1 = 1.0f / d1; // Use reciprocal to avoid repeated division
        double term1 = (att1 - 1.0f) * inv_d1;
        
        intensity_mu1 = intensity_mu1 * att1 + (1.0f + term1) * p_next + (-att1 - term1) * p_curr;

        // Calculation for Mu 2
        double d2 = tau_layer / mu2;
        double att2 = exp(-d2);
        double inv_d2 = 1.0f / d2;
        double term2 = (att2 - 1.0f) * inv_d2;
                
        intensity_mu2 = intensity_mu2 * att2 + (1.0f + term2) * p_next + (-att2 - term2) * p_curr;
      }

      // Propagate values to next iteration without re-reading/re-calculating
      z_curr = z_next;
      p_curr = p_next;
    }

    // Final weighted integration
    // Combined constants: 2 * pi * 0.5 * 1e-3 = pi * 1e-3
    const double final_const = constants::pi * 1e-3 * spectrum_scaling;
    model_spectrum_gpu[tid] = final_const * (intensity_mu1 * mu1 + intensity_mu2 * mu2);
  }
}



__global__ 
void shortCharacteristicsDev_Shared(
  double* __restrict__ model_spectrum_gpu,
  const float* __restrict__ absorption_coeff_dev,
  const double* __restrict__ wavenumber_list_dev,
  const float* __restrict__ cloud_optical_depth_dev,
  const float* __restrict__ temperature_dev,
  const float* __restrict__ vertical_grid_dev,
  const double spectrum_scaling,
  const int nb_spectral_points,
  const int nb_grid_points)
{
  // Dynamically or statically allocate shared memory
  // Note: Adjust the size or use extern __shared__ if nb_grid_points is large
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
    
    for (int i = 0; i < nb_grid_points - 1; ++i)
    {
      const float z_next = s_vertical_grid[i+1];
      const float p_next = planckFunction(s_temperature[i+1], wn_cube, wavenumber);
      
      // Global memory loads (absorption is unique per spectral point/thread)
      float abs_curr = absorption_coeff_dev[i * nb_spectral_points + tid];
      float abs_next = absorption_coeff_dev[(i + 1) * nb_spectral_points + tid];
      
      float tau_layer = (z_next - z_curr) * (abs_next + abs_curr) * 0.5;
      
      if (cloud_optical_depth_dev != nullptr) 
        tau_layer += cloud_optical_depth_dev[i * nb_spectral_points + tid];
      
      if (tau_layer > 1e-12) 
      {
        // Optimization: Pre-calculate shared terms
        double inv_mu1 = 1.0f / mu1;
        double inv_mu2 = 1.0f / mu2;
        
        // Mu 1 Path
        double d1 = tau_layer * inv_mu1;
        double att1 = exp(-d1);
        double term1 = (att1 - 1.0f) / d1;
        intensity_mu1 = intensity_mu1 * att1 + (1.0f + term1) * p_next + (-att1 - term1) * p_curr;
        
        // Mu 2 Path
        double d2 = tau_layer * inv_mu2;
        double att2 = exp(-d2);
        double term2 = (att2 - 1.0f) / d2;
        intensity_mu2 = intensity_mu2 * att2 + (1.0f + term2) * p_next + (-att2 - term2) * p_curr;
      }
      
      z_curr = z_next;
      p_curr = p_next;
    }

    const double final_const = 3.141592653589793 * 1e-3 * spectrum_scaling;
    model_spectrum_gpu[tid] = final_const * (intensity_mu1 * mu1 + intensity_mu2 * mu2);
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
  double* model_spectrum_dev)
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

  // shortCharacteristicsDev<<<blocks,threads>>>(
  //   model_spectrum_dev,
  //   absorption_coeff_dev,
  //   spectral_grid->wavenumber_list_gpu,
  //   cloud_optical_depth_dev,
  //   atmosphere.temperature_dev,
  //   atmosphere.altitude_dev,
  //   spectrum_scaling,
  //   nb_spectral_points,
  //   nb_grid_points);

  // shortCharacteristicsDevOld<<<blocks,threads>>>(
  //   model_spectrum_dev,
  //   absorption_coeff_dev,
  //   spectral_grid->wavenumber_list_gpu,
  //   cloud_optical_depth_dev,
  //   atmosphere.temperature_dev,
  //   atmosphere.altitude_dev,
  //   spectrum_scaling,
  //   nb_spectral_points,
  //   nb_grid_points);


  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
}


}
