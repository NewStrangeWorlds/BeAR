/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2020 Daniel Kitzmann
*
* Helios-r2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Helios-r2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* Helios-r2 directory under <LICENSE>. If not, see
* <http://www.gnu.org/licenses/>.
*/


#include "short_characteristics_kernels.h"


#include <iostream>
#include <vector>
#include "math.h"
#include <stdio.h>
#include <new>

#include "data_management_kernels.h"
#include "../additional/physical_const.h"


#include "error_check.h"
#include "planck_function.h"

namespace helios{


//solves the radiative transfer equation with the short characteristic method
//uses two angles, distributed according to a Gaussian quadrature scheme
__global__ void shortCharacteristicsDevice(double* model_spectrum_gpu,
                                           const double* absorption_coeff_dev, const double* wavenumber_list_dev,
                                           const double* temperature_dev, const double* vertical_grid_dev,
                                           const double radius_distance_scaling,
                                           const int nb_spectral_points, const int nb_grid_points)
{
  //tid is the wavenumber index
  //each core calculates one wavenumber
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_spectral_points; tid += blockDim.x * gridDim.x)
  {
    //inner boundary condition
    double intensity_mu1 = planckFunction(temperature_dev[0], wavenumber_list_dev[tid]);
    double intensity_mu2 = intensity_mu1;


    for (int i=0; i<nb_grid_points-1; ++i)
    {
      const double planck_function1 = planckFunction(temperature_dev[i], wavenumber_list_dev[tid]);
      const double planck_function2 = planckFunction(temperature_dev[i+1], wavenumber_list_dev[tid]);
      const double optical_depth_layer = (vertical_grid_dev[i+1] - vertical_grid_dev[i]) 
                                        * ( absorption_coeff_dev[(i+1)*nb_spectral_points + tid ] + absorption_coeff_dev[i*nb_spectral_points + tid])/2.;
      
      //Gauss angle 1
      const double delta1 = optical_depth_layer/0.339981;
      const double attenuation_factor1 = exp(-delta1);

      intensity_mu1 = intensity_mu1 * attenuation_factor1;

      const double beta1 = 1.0 + (attenuation_factor1 - 1.0)/delta1;
      const double gamma1 = -attenuation_factor1 - (attenuation_factor1 - 1.0)/delta1;

      intensity_mu1 += beta1 * planck_function2 + gamma1 * planck_function1;


      //Gauss angle 2
      const double delta2 = optical_depth_layer/0.861136;
      const double attenuation_factor2 = exp(-delta2);

      intensity_mu2 = intensity_mu2 * attenuation_factor2;

      const double beta2 = 1.0 + (attenuation_factor2 - 1.0)/delta2;
      const double gamma2 = -attenuation_factor2 - (attenuation_factor2 - 1.0)/delta2;

      intensity_mu2 += beta2 * planck_function2 + gamma2 * planck_function1;
    }

    //and integration with the corresponding Gauss-Legendre quadrature weights
    model_spectrum_gpu[tid] = 2.0 * constants::pi * (intensity_mu1 * 0.339981 * 0.652145 + intensity_mu2 * 0.861136 * 0.347855) * 1e-3 * radius_distance_scaling; //in W m-2 cm-1
  }

}



//solves the radiative transfer equation with the short characteristic method
//uses two angles, distributed according to a Gaussian quadrature scheme
//additionally includes a cloud layer
__global__ void shortCharacteristicsCloudDevice(double* model_spectrum_gpu,
                                                const double* absorption_coeff_dev, const double* wavenumber_list_dev,
                                                const double* cloud_optical_depth_dev,
                                                const double* temperature_dev, const double* vertical_grid_dev,
                                                const double radius_distance_scaling,
                                                const int nb_spectral_points, const int nb_grid_points)
{
  //tid is the wavenumber index
  //int tid = blockIdx.x * blockDim.x + threadIdx.x;
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_spectral_points; tid += blockDim.x * gridDim.x)
  {
    //inner boundary condition
    double intensity_mu1 = planckFunction(temperature_dev[0], wavenumber_list_dev[tid]);
    double intensity_mu2 = intensity_mu1;


    for (int i=0; i<nb_grid_points-1; ++i)
    {
      const double planck_function1 = planckFunction(temperature_dev[i], wavenumber_list_dev[tid]);
      const double planck_function2 = planckFunction(temperature_dev[i+1], wavenumber_list_dev[tid]);

      //optical depth includes the molecular absorption coefficients as well as the one of the cloud
      const double optical_depth_layer = (vertical_grid_dev[i+1] - vertical_grid_dev[i]) 
                                       * ( absorption_coeff_dev[(i+1)*nb_spectral_points + tid ] + absorption_coeff_dev[i*nb_spectral_points + tid])/2.
                                         + cloud_optical_depth_dev[i];

      //Gauss angle 1
      const double delta1 = optical_depth_layer/0.339981;
      const double attenuation_factor1 = exp(-delta1);

      intensity_mu1 = intensity_mu1 * attenuation_factor1;
 
      const double beta1 = 1.0 + (attenuation_factor1 - 1.0)/delta1;
      const double gamma1 = -attenuation_factor1 - (attenuation_factor1 - 1.0)/delta1;

      intensity_mu1 += beta1 * planck_function2 + gamma1 * planck_function1;


      //Gauss angle 2
      const double delta2 = optical_depth_layer/0.861136;
      const double attenuation_factor2 = exp(-delta2);

      intensity_mu2 = intensity_mu2 * attenuation_factor2;

      const double beta2 = 1.0 + (attenuation_factor2 - 1.0)/delta2;
      const double gamma2 = -attenuation_factor2 - (attenuation_factor2 - 1.0)/delta2;

      intensity_mu2 += beta2 * planck_function2 + gamma2 * planck_function1;
    }
    
    //and integration with the corresponding Gauss-Legendre quadrature weights
    model_spectrum_gpu[tid] = 2.0 * constants::pi * (intensity_mu1 * 0.339981 * 0.652145 + intensity_mu2 * 0.861136 * 0.347855) * 1e-3 * radius_distance_scaling; //in W m-2 cm-1
  }

}





//This the version including an optional cloud layer
__host__ void shortCharacteristicsGPU(double* model_spectrum_dev,
                                      double* absorption_coeff_dev, double* wavenumber_list_dev,
                                      const std::vector<double>& cloud_optical_depth,
                                      const std::vector<double>& temperature, const std::vector<double>& vertical_grid,
                                      const double radius_distance_scaling,
                                      const size_t& nb_spectral_points)
{
  size_t nb_grid_points = temperature.size();

  const int bytes = nb_grid_points*sizeof(double);



  double* temperature_dev = nullptr;

  cudaMalloc(&temperature_dev, bytes);
  cudaMemcpy(temperature_dev, &temperature[0], bytes, cudaMemcpyHostToDevice);


  double* vertical_grid_dev = nullptr;

  cudaMalloc(&vertical_grid_dev, bytes);
  cudaMemcpy(vertical_grid_dev, &vertical_grid[0], bytes, cudaMemcpyHostToDevice);


  double* cloud_optical_depth_dev = nullptr;

  cudaMalloc(&cloud_optical_depth_dev, bytes);
  cudaMemcpy(cloud_optical_depth_dev, &cloud_optical_depth[0], bytes-1, cudaMemcpyHostToDevice);


  cudaThreadSynchronize();


  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;



  shortCharacteristicsCloudDevice<<<blocks,threads>>>(model_spectrum_dev,
                                                      absorption_coeff_dev, wavenumber_list_dev,
                                                      cloud_optical_depth_dev,
                                                      temperature_dev, vertical_grid_dev,
                                                      radius_distance_scaling,
                                                      nb_spectral_points, nb_grid_points);


  cudaThreadSynchronize();

  cudaFree(temperature_dev);
  cudaFree(vertical_grid_dev);
  cudaFree(cloud_optical_depth_dev);

  gpuErrchk( cudaPeekAtLastError() );
}



__host__ void shortCharacteristicsGPU(double* model_spectrum_dev,
                                      double* absorption_coeff_dev, double* wavenumber_list_dev,
                                      const std::vector<double>& temperature, const std::vector<double>& vertical_grid,
                                      const double radius_distance_scaling,
                                      const size_t& nb_spectral_points)
{
  size_t nb_grid_points = temperature.size();

  const int bytes = nb_grid_points*sizeof(double);



  double* temperature_dev = nullptr;

  cudaMalloc(&temperature_dev, bytes);
  cudaMemcpy(temperature_dev, &temperature[0], bytes, cudaMemcpyHostToDevice);


  double* vertical_grid_dev = nullptr;

  cudaMalloc(&vertical_grid_dev, bytes);
  cudaMemcpy(vertical_grid_dev, &vertical_grid[0], bytes, cudaMemcpyHostToDevice);


  cudaThreadSynchronize();


  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;



  shortCharacteristicsDevice<<<blocks,threads>>>(model_spectrum_dev,
                                                 absorption_coeff_dev, wavenumber_list_dev,
                                                 temperature_dev, vertical_grid_dev,
                                                 radius_distance_scaling,
                                                 nb_spectral_points, nb_grid_points);


  cudaThreadSynchronize();

  cudaFree(temperature_dev);
  cudaFree(vertical_grid_dev);

  gpuErrchk( cudaPeekAtLastError() );
}





}
