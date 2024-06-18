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


namespace helios{

//solves the radiative transfer equation with the short characteristic method
//uses two angles, distributed according to a Gaussian quadrature scheme
__global__ void shortCharacteristicsDev(
  double* model_spectrum_gpu,
  const double* absorption_coeff_dev,
  const double* wavenumber_list_dev,
  const double* cloud_optical_depth_dev,
  const double* temperature_dev,
  const double* vertical_grid_dev,
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
        (vertical_grid_dev[i+1] - vertical_grid_dev[i]) 
        * ( absorption_coeff_dev[(i+1)*nb_spectral_points + tid ] + absorption_coeff_dev[i*nb_spectral_points + tid])/2.;

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



void ShortCharacteristics::calcSpectrumGPU(
  const Atmosphere& atmosphere,
  double* absorption_coeff_dev,
  double* scattering_coeff_dev,
  double* cloud_optical_depth_dev,
  double* cloud_single_scattering_dev,
  double* cloud_asym_param_dev,
  const double spectrum_scaling,
  double* model_spectrum_dev)
{
  size_t nb_grid_points = atmosphere.temperature.size();
  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();


  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;


  shortCharacteristicsDev<<<blocks,threads>>>(
    model_spectrum_dev,
    absorption_coeff_dev,
    spectral_grid->wavenumber_list_gpu,
    cloud_optical_depth_dev,
    atmosphere.temperature_dev,
    atmosphere.altitude_dev,
    spectrum_scaling,
    nb_spectral_points,
    nb_grid_points);


  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
}


}
