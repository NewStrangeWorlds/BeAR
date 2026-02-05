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


#include "../forward_model/transmission/transmission.h"
#include "../forward_model/atmosphere/atmosphere.h"

#include <iostream>
#include <cmath>
#include <stdio.h>

#include "error_check.h"
#include "reduce_kernels.h"
#include "../additional/physical_const.h"


namespace bear{


__global__
void transmissionSpectrumKernel(
  const int nb_spectral_points,
  const int nb_grid_points,
  const double bottom_radius,
  const double star_radius,
  const double* __restrict__ altitude,
  const float* __restrict__ extinction_coeff,
  double* spectrum)
{
  //thread index is the wavelength index
  int w = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (w >= nb_spectral_points) return;

  //each thread integrates a single wavelength
  float effective_tangent_height = 0.0;

  //integrate effective tangent height
  float prev_transmission = 1.0;  // top limb = transparent
  float prev_alt   = altitude[nb_grid_points-1];

  for (int t = nb_grid_points-2; t >= 0; --t)
  {
    //tangent optical depth
    float tau = 0.0;
    float b = bottom_radius + altitude[t];

    for (int i = t; i < nb_grid_points-1; ++i)
    {
      float a1 = bottom_radius + altitude[i];
      float a2 = bottom_radius + altitude[i+1];
      
      float path = sqrt(a2*a2 - b*b);
      
      if (i != t)
        path -= sqrt(a1*a1-b*b);
      
      int idx1 = i*nb_spectral_points + w;
      int idx2 = (i+1)*nb_spectral_points + w;

      float ext1 = extinction_coeff[idx1];
      float ext2 = extinction_coeff[idx2];
      
      tau += path * (ext1 + ext2);
      
      if (tau > transmission_optical_depth_cutoff) break;
    }
    
    float transmission = exp(-tau);

    //trapezoidal integration for effective tangent height
    //the factor of 0.5 from the trapezoidal rule cancels with 
    //the factor of 2 from accounting for both hemispheres
    float alt = altitude[t];

    effective_tangent_height +=
      ( (bottom_radius + alt) * (1.0 - transmission)
      + (bottom_radius + prev_alt) * (1.0 - prev_transmission) )
      * (prev_alt - alt);

    prev_transmission = transmission;
    prev_alt   = alt;
    
    //if the optical depth is too high, we simply add the bottom part
    if (transmission < transmission_cutoff)
    {
      //effective_tangent_height += (bottom_radius + alt)*(bottom_radius + alt) - bottom_radius*bottom_radius;
      effective_tangent_height += (2*bottom_radius + alt)*alt;
      break;
    }
  }
  
  effective_tangent_height =
      sqrt(effective_tangent_height + bottom_radius*bottom_radius)
      - bottom_radius;

  double planet_radius = effective_tangent_height + bottom_radius;

  spectrum[w] = (planet_radius * planet_radius) /
                (star_radius * star_radius) * 1e6;
}



//sum all transparencies to get the total extinction coefficient
//the result is placed in the absorption_coeff array
__global__ 
void sumExtinctionCoeff(
  const int nb_points, 
  float* __restrict__ absorption_coeff, 
  const float* __restrict__ scattering_coeff)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;  // Global thread index
  
  if (i < nb_points)
    absorption_coeff[i] += scattering_coeff[i];
}


//sum all transparencies to get the total extinction coefficient
//the result is placed in the absorption_coeff array
__global__ 
void sumExtinctionCoeff(
  const int nb_points, 
  float* __restrict__ absorption_coeff, 
  const float* __restrict__ scattering_coeff,
  const float* __restrict__ cloud_extinction_coeff)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;  // Global thread index
  
  if (i < nb_points)
    absorption_coeff[i] += scattering_coeff[i] + cloud_extinction_coeff[i];
}



__host__ 
void  TransmissionModel::calcTransitDepthGPU(
  double* transit_radius_dev, 
  float* absorption_coeff_dev, 
  float* scattering_coeff_dev, 
  float* cloud_extinction_dev, 
  const Atmosphere& atmosphere, 
  const size_t nb_spectral_points, 
  const double radius_planet, 
  const double radius_star)
{
  int threads = 256;
  int blocks1 = (nb_spectral_points*nb_grid_points) / threads;

  if ((nb_spectral_points*nb_grid_points) % threads) blocks1++;
  
  if (cloud_extinction_dev != nullptr)
    sumExtinctionCoeff<<<blocks1,threads>>>(
      nb_spectral_points*nb_grid_points, 
      absorption_coeff_dev, 
      scattering_coeff_dev, 
      cloud_extinction_dev); 
  else
    sumExtinctionCoeff<<<blocks1,threads>>>(
      nb_spectral_points*nb_grid_points, 
      absorption_coeff_dev, 
      scattering_coeff_dev);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );


  int blocks  = (nb_spectral_points + threads - 1) / threads;

  transmissionSpectrumKernel<<<blocks, threads>>>(
    nb_spectral_points,
    nb_grid_points,
    radius_planet,
    radius_star,
    atmosphere.altitude_dev,
    absorption_coeff_dev,
    transit_radius_dev);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}

}