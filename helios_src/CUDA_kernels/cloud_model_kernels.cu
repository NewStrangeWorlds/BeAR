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


#include "../cloud_model/grey_cloud_model.h"
#include "../cloud_model/kh_cloud_model.h"
#include "../spectral_grid/spectral_grid.h"

#include <iostream>
#include <cmath>
#include <stdio.h>

#include "error_check.h"
#include "reduce_kernels.h"
#include "../additional/physical_const.h"


namespace helios{


__global__ void greyCloudModel(
  double* optical_depth, 
  double* single_scattering_albedo, 
  const double layer_optical_depth,
  const int cloud_top,
  const int cloud_bottom,
  const int nb_spectral_points)
{
  //the thread index tid is the wavelength
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_spectral_points; tid += blockDim.x * gridDim.x)
  {
    for (int i=cloud_top; i>cloud_bottom; --i)
    { 
      const int j = i*nb_spectral_points + tid;
      const double optical_depth_mix = optical_depth[j] + layer_optical_depth;

      single_scattering_albedo[j] = optical_depth[j] * single_scattering_albedo[j] / optical_depth_mix;

      optical_depth[j] = optical_depth_mix;
    }
  }

}



__global__ void KHnongreyCloudModel(
  double* optical_depth, 
  double* single_scattering_albedo, 
  double* wavelengths,
  const double optical_depth_layer_ref,
  const double reference_q,
  const double particle_size,
  const double a0,
  const double q0,
  const int cloud_top,
  const int cloud_bottom,
  const int nb_spectral_points)
{
  //the thread index tid is the wavelength
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_spectral_points; tid += blockDim.x * gridDim.x)
  {
    const double size_parameter = 2*constants::pi * particle_size / wavelengths[tid];

    const double optical_depth_layer = optical_depth_layer_ref * reference_q 
      / (q0 * pow(size_parameter, -a0) + pow(size_parameter, 0.2) );

    for (int i=cloud_top; i>cloud_bottom; --i)
    {
      const int j = i*nb_spectral_points + tid;
      const double optical_depth_mix = optical_depth[j] + optical_depth_layer;

      single_scattering_albedo[j] = optical_depth[j] * single_scattering_albedo[j] / optical_depth_mix;
      
      optical_depth[j] = optical_depth_mix;
    }
  }

}




//calculates the vertical distribution of the grey layer
//needs three parameters: cloud top pressure, cloud bottom (fraction of top pressure), and optical depth
//the optical depth will be distributed over the layers between the cloud's top and bottom
__host__ void GreyCloudModel::opticalPropertiesGPU(const std::vector<double>& parameters, const Atmosphere& atmosphere,
  SpectralGrid* spectral_grid,
  double* optical_depth_dev, 
  double* single_scattering_dev, 
  double* asym_param)
{
  double cloud_optical_depth = parameters[0];
  double cloud_top_pressure = parameters[1];
  double cloud_bottom_pressure = 0;

  if (fixed_bottom == false)
    cloud_bottom_pressure = cloud_top_pressure * parameters[2];


  if (cloud_optical_depth < 0) cloud_optical_depth = 0;


  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();
  size_t nb_grid_points = atmosphere.nb_grid_points;

  unsigned int cloud_top_index = 0;
  unsigned int cloud_bottom_index = 0;

  if (fixed_bottom == true)
    cloudPosition(atmosphere, cloud_top_pressure, cloud_top_index, cloud_bottom_index);
  else
    cloudPosition(atmosphere, cloud_top_pressure, cloud_bottom_pressure, cloud_top_index, cloud_bottom_index);


  const double optical_depth_layer = cloud_optical_depth/static_cast<double>(cloud_top_index - cloud_bottom_index);

  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  greyCloudModel<<<blocks,threads>>>(
    optical_depth_dev,
    single_scattering_dev,
    optical_depth_layer,
    cloud_top_index,
    cloud_bottom_index,
    nb_spectral_points);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}



//calculates the vertical distribution of the grey layer
//needs three parameters: cloud top pressure, cloud bottom (fraction of top pressure), and optical depth
//the optical depth will be distributed over the layers between the cloud's top and bottom
__host__ void KHCloudModel::opticalPropertiesGPU(const std::vector<double>& parameters, const Atmosphere& atmosphere,
  SpectralGrid* spectral_grid,
  double* optical_depth_dev, 
  double* single_scattering_dev, 
  double* asym_param)
{
  double cloud_optical_depth = parameters[0];
  double q0 = parameters[1];
  double a0 = parameters[2];
  double particle_size = parameters[3];
  double cloud_top_pressure = parameters[4];


  double cloud_bottom_pressure = 0;

  if (fixed_bottom == false)
    cloud_bottom_pressure = cloud_top_pressure * parameters[5];


  if (cloud_optical_depth < 0) cloud_optical_depth = 0;


  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();
  size_t nb_grid_points = atmosphere.nb_grid_points;

  unsigned int cloud_top_index = 0;
  unsigned int cloud_bottom_index = 0;

  if (fixed_bottom == true)
    cloudPosition(atmosphere, cloud_top_pressure, cloud_top_index, cloud_bottom_index);
  else
    cloudPosition(atmosphere, cloud_top_pressure, cloud_bottom_pressure, cloud_top_index, cloud_bottom_index);


  double optical_depth_layer_reference = cloud_optical_depth/static_cast<double>(cloud_top_index - cloud_bottom_index);

  double reference_size_parameter = 2*constants::pi * particle_size / reference_wavelength;
  double reference_q = q0 * std::pow(reference_size_parameter, -a0) + std::pow(reference_size_parameter, 0.2);

  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  KHnongreyCloudModel<<<blocks,threads>>>(
    optical_depth_dev,
    single_scattering_dev,
    spectral_grid->wavelength_list_gpu,
    optical_depth_layer_reference,
    reference_q,
    particle_size,
    a0,
    q0,
    cloud_top_index,
    cloud_bottom_index,
    nb_spectral_points);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}


}
