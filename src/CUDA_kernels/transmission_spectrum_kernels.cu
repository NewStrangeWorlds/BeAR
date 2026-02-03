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

__forceinline__ __device__ double distanceToTangentCenter(
  const double tangent_altitude, 
  const double altitude, 
  const double bottom_radius)
{
  const double a = bottom_radius + altitude;
  const double b = bottom_radius + tangent_altitude;

  if (a <= b) return 0;

  return sqrt(a*a - b*b);
}

__global__
void transmissionSpectrumKernel(
  int nb_spectral_points,
  int nb_grid_points,
  double bottom_radius,
  double star_radius,
  const double* __restrict__ altitude,
  const double* __restrict__ absorption,
  const double* __restrict__ scattering,
  const double* __restrict__ cloud,
  double transmission_cutoff,
  double* spectrum)
{
  int w = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (w >= nb_spectral_points) return;

  // Each thread gets its own local transmission array (in registers/local memory)
  double effective_tangent_height = 0.0;

  // ---- integrateEffectiveTangentHeight ----
  double prev_trans = 1.0;  // top limb = transparent
  double prev_alt   = altitude[nb_grid_points-1];

  for (int t = nb_grid_points-2; t >= 0; --t)
  {
    // ---- tangentOpticalDepth ----
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

      // int idx1 = w*nb_grid_points + i;
      // int idx2 = w*nb_grid_points + i+1;

      float ext1 = absorption[idx1] + scattering[idx1]; // + cloud[idx1];
      float ext2 = absorption[idx2] + scattering[idx2]; // + cloud[idx2];

      tau += path * (ext1 + ext2);
      
      if (tau > transmission_cutoff) break;
    }
    
    float trans = exp(-tau);

    // ---- trapezoidal integration for effective tangent height ----
    double alt = altitude[t];

    effective_tangent_height +=
      ( (bottom_radius + alt) * (1.0 - trans)
      + (bottom_radius + prev_alt) * (1.0 - prev_trans) )
      * (prev_alt - alt);

    prev_trans = trans;
    prev_alt   = alt;
  }
  
  effective_tangent_height =
      sqrt(effective_tangent_height + bottom_radius*bottom_radius)
      - bottom_radius;

  double planet_radius = effective_tangent_height + bottom_radius;

  spectrum[w] = (planet_radius * planet_radius) /
                (star_radius * star_radius) * 1e6;
  
}



__forceinline__ __device__ double tangentOpticalDepth(
  const double tangent_radius, 
  const double altitude1, 
  const double altitude2, 
  const double extinction_coeff1, 
  const double extinction_coeff2, 
  const double bottom_radius)
{
  const double path_length = distanceToTangentCenter(tangent_radius, altitude2, bottom_radius) 
                           - distanceToTangentCenter(tangent_radius, altitude1, bottom_radius);

  //account for both hemispheres by multiplying the results by 2
  //note: this cancels the factor of 0.5 from the trapezoidal rule
  return path_length * (extinction_coeff1 + extinction_coeff2);
}



__forceinline__ __device__ double tangentPathsTransmission(
  const int wavelength_index, 
  const int tangent_radius_index, 
  const double bottom_radius, 
  const int nb_grid_points, 
  const int nb_spectral_points, 
  double* altitudes, 
  double* extinction_coeff)
{
  double optical_depth = 0;
  double altitude1 = altitudes[tangent_radius_index]+0.01;
  double altitude2 = 0;
  const double tangent_radius = altitudes[tangent_radius_index];

  double extinction_coeff1 = extinction_coeff[(tangent_radius_index)*nb_spectral_points + wavelength_index];
  double extinction_coeff2;

  for (int i=tangent_radius_index; i<nb_grid_points-1; ++i)
  {
    altitude2 = altitudes[i+1];

    extinction_coeff2 = extinction_coeff[(i+1)*nb_spectral_points + wavelength_index];

    optical_depth += tangentOpticalDepth(
      tangent_radius, 
      altitude1, 
      altitude2, 
      extinction_coeff1, 
      extinction_coeff2, 
      bottom_radius);

    altitude1 = altitude2;
    extinction_coeff1 = extinction_coeff2;

    if (optical_depth > transmission_optical_depth_cutoff) break;
  }

  return exp(-optical_depth);
}



__device__ double effectiveTangentHeight(
  const int wavelength_index, 
  const int nb_grid_points, 
  const int nb_spectral_points, 
  double* altitudes, 
  double* extinction_coeff,
  const double radius_planet)
{
  double effective_tangent_height = 0;
  double path_transmission1 = 1;
  double path_transmission2 = 0;

  double altitude1 = altitudes[nb_grid_points-1];
  double altitude2 = 0;

  for (int i=nb_grid_points-1; i>0; --i)
  {
    altitude2 = altitudes[i-1];

    path_transmission2 = tangentPathsTransmission(
      wavelength_index, 
      i-1, 
      radius_planet, 
      nb_grid_points, 
      nb_spectral_points, 
      altitudes, 
      extinction_coeff);

    effective_tangent_height += 2. * ( (radius_planet + altitude1) * (1. - path_transmission1)
                                     + (radius_planet + altitude2) * (1. - path_transmission2) )
                                   * (altitude1 - altitude2) * 0.5;

    if (path_transmission2 < transmission_cutoff)
    {
      effective_tangent_height += ((radius_planet + altitude2) + radius_planet)*altitude2;
      break;
    }

    altitude1 = altitude2;
    path_transmission1 = path_transmission2;
  }


  effective_tangent_height = sqrt(effective_tangent_height + radius_planet*radius_planet) - radius_planet;
 
  return effective_tangent_height;
}


//sum all transparencies to get the total extinction coefficient
//the result is placed in the absorption_coeff array
__global__ void sumExtinctionCoeff(
  const int nb_grid_points, 
  const int nb_spectral_points, 
  double* absorption_coeff, 
  double* scattering_coeff,
  double* cloud_extinction_coeff)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < nb_spectral_points; i += blockDim.x * gridDim.x)
  {
    for (int j = 0; j < nb_grid_points; ++j)
    {
      absorption_coeff[j*nb_spectral_points + i] += scattering_coeff[j*nb_spectral_points + i] + cloud_extinction_coeff[j*nb_spectral_points + i];
    }
  }
}


//sum all transparencies to get the total extinction coefficient
//the result is placed in the absorption_coeff array
__global__ void sumExtinctionCoeff(
  const int nb_grid_points, 
  const int nb_spectral_points, 
  double* absorption_coeff, 
  double* scattering_coeff)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < nb_spectral_points; i += blockDim.x * gridDim.x)
  {
    for (int j = 0; j < nb_grid_points; ++j)
    {
      absorption_coeff[j*nb_spectral_points + i] += scattering_coeff[j*nb_spectral_points + i];
    }
  }
}


//sum all transparencies to get the total extinction coefficient
//the result is placed in the absorption_coeff array
__global__ void sumExtinctionCoeff(
  const int nb_points, 
  double* absorption_coeff, 
  double* scattering_coeff)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;  // Global thread index
  
  if (i < nb_points)
    absorption_coeff[i] += scattering_coeff[i];
}


__global__ void transitDepthDev(
  double* transit_radius,
  double* extinction_coeff,
  double* altitudes, 
  const int nb_spectral_points, 
  const int nb_grid_points,
  const double radius_planet, 
  const double radius_star_sqr)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < nb_spectral_points; i += blockDim.x * gridDim.x)
  {
    transit_radius[i] = effectiveTangentHeight(
      i, 
      nb_grid_points, 
      nb_spectral_points, 
      altitudes, 
      extinction_coeff,
      radius_planet) + radius_planet;

    transit_radius[i] = transit_radius[i]*transit_radius[i]/(radius_star_sqr) * 1e6;
  }

}



__host__ void  TransmissionModel::calcTransitDepthGPU(
  double* transit_radius_dev, 
  double* absorption_coeff_dev, 
  double* scattering_coeff_dev, 
  double* cloud_extinction_dev, 
  const Atmosphere& atmosphere, 
  const size_t nb_spectral_points, 
  const double radius_planet, 
  const double radius_star)
{
  // int threads = 256;
  // int blocks = (nb_spectral_points*nb_grid_points) / threads;

  // if ((nb_spectral_points*nb_grid_points) % threads) blocks++;

  
  // if (cloud_extinction_dev != nullptr)
  //   sumExtinctionCoeff<<<blocks,threads>>>(
  //     nb_spectral_points, 
  //     atmosphere.nb_grid_points, 
  //     absorption_coeff_dev, 
  //     scattering_coeff_dev, 
  //     cloud_extinction_dev); 
  // else
  //   // sumExtinctionCoeff<<<blocks,threads>>>(
  //   //   nb_spectral_points, 
  //   //   atmosphere.nb_grid_points, 
  //   //   absorption_coeff_dev, 
  //   //   scattering_coeff_dev);
  //    sumExtinctionCoeff<<<blocks,threads>>>(
  //     nb_spectral_points*nb_grid_points, 
  //     absorption_coeff_dev, 
  //     scattering_coeff_dev);

  // cudaDeviceSynchronize();
  // gpuErrchk( cudaPeekAtLastError() );
  // gpuErrchk( cudaDeviceSynchronize() );


  // threads = 256;
  // blocks = nb_spectral_points / threads;

  // if (nb_spectral_points % threads) blocks++;

  // transitDepthDev<<<blocks,threads>>>(
  //   transit_radius_dev, 
  //   absorption_coeff_dev, 
  //   atmosphere.altitude_dev, 
  //   nb_spectral_points, 
  //   nb_grid_points, 
  //   radius_planet, 
  //   radius_star*radius_star);


  // cudaDeviceSynchronize();
  // gpuErrchk( cudaPeekAtLastError() );
  // gpuErrchk( cudaDeviceSynchronize() );

  int threads = 256;
  int blocks  = (nb_spectral_points + threads - 1) / threads;

  transmissionSpectrumKernel<<<blocks, threads>>>(
    nb_spectral_points,
    nb_grid_points,
    radius_planet,
    radius_star,
    atmosphere.altitude_dev,
    absorption_coeff_dev,
    scattering_coeff_dev,
    cloud_extinction_dev,
    transmission_optical_depth_cutoff,
    transit_radius_dev);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}

}