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


#include "../forward_model/transmission/transmission.h"
#include "../forward_model/atmosphere/atmosphere.h"

#include <iostream>
#include <cmath>
#include <stdio.h>

#include "error_check.h"
#include "reduce_kernels.h"
#include "../additional/physical_const.h"


namespace helios{


__forceinline__ __device__ double distanceToTangentDev(const double tangent_altitude, const double altitude, const double bottom_radius)
{
  const double a = bottom_radius + altitude;
  const double b = bottom_radius + tangent_altitude;

  return sqrt(a*a - b*b);
}



__forceinline__ __device__ double tangentOpticalDepthDev(
  const double tangent_radius, 
  const double altitude1, 
  const double altitude2, 
  const double extinction_coeff1, 
  const double extinction_coeff2, 
  const double bottom_radius)
{
  const double path_length = distanceToTangentDev(tangent_radius, altitude2, bottom_radius) -  distanceToTangentDev(tangent_radius, altitude1, bottom_radius);
  
  return path_length * (extinction_coeff1 + extinction_coeff2);
}



__forceinline__ __device__ double tangentPathsTransmissionDev(
  const int wavelength_index, 
  const int tangent_radius_index, 
  const double bottom_radius, 
  const int nb_grid_points, 
  const int nb_spectral_points, 
  double* altitudes, 
  double* absorption_coeff_dev, 
  double* scattering_coeff_dev)
{
  double optical_depth = 0;
  double altitude1 = altitudes[tangent_radius_index]+0.01;
  double altitude2 = 0;
  const double tangent_radius = altitudes[tangent_radius_index];

  double extinction_coeff1 = absorption_coeff_dev[(tangent_radius_index)*nb_spectral_points + wavelength_index] + scattering_coeff_dev[(tangent_radius_index)*nb_spectral_points + wavelength_index];
  double extinction_coeff2;

  for (int i=tangent_radius_index; i<nb_grid_points-1; ++i)
  {
    altitude2 = altitudes[i+1];

    extinction_coeff2 = absorption_coeff_dev[(i+1)*nb_spectral_points + wavelength_index] + scattering_coeff_dev[(i+1)*nb_spectral_points + wavelength_index];

    optical_depth += tangentOpticalDepthDev(tangent_radius, altitude1, altitude2, extinction_coeff1, extinction_coeff2, bottom_radius);

    altitude1 = altitude2;
    extinction_coeff1 = extinction_coeff2;


    if (optical_depth > transmission_optical_depth_cutoff) break;
  }

  return exp(-optical_depth);
}



__device__ double effectiveTangentHeightDev(
  const int wavelength_index, 
  const int nb_grid_points, 
  const int nb_spectral_points, 
  double* altitudes, 
  double* absorption_coeff_dev, 
  double* scattering_coeff_dev, 
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

    path_transmission2 = tangentPathsTransmissionDev(
      wavelength_index, 
      i-1, 
      radius_planet, 
      nb_grid_points, 
      nb_spectral_points, 
      altitudes, 
      absorption_coeff_dev, 
      scattering_coeff_dev);

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



__global__ void transitDepthDevice(
  double* transit_radius_dev,
  double* absorption_coeff_dev, 
  double* scattering_coeff_dev,
  double* altitudes, 
  const int nb_spectral_points, 
  const int nb_grid_points,
  const double radius_planet, 
  const double radius_star_sqr)
{
  
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < nb_spectral_points; i += blockDim.x * gridDim.x)
  {

   transit_radius_dev[i] = effectiveTangentHeightDev(i, nb_grid_points, nb_spectral_points, altitudes, absorption_coeff_dev, scattering_coeff_dev, radius_planet) + radius_planet;
   transit_radius_dev[i] = transit_radius_dev[i]*transit_radius_dev[i]/(radius_star_sqr) * 1e6;
  }

}



__host__ void  TransmissionModel::calcTransitRadiusGPU(
  double* transit_radius_dev, 
  double* absorption_coeff_dev, 
  double* scattering_coeff_dev, 
  const Atmosphere& atmosphere, 
  const size_t nb_spectral_points, 
  const double radius_planet, 
  const double radius_star)
{

  int threads = 256;

  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;


  transitDepthDevice<<<blocks,threads>>>(transit_radius_dev, absorption_coeff_dev, scattering_coeff_dev, atmosphere.altitude_dev, nb_spectral_points, nb_grid_points, radius_planet, radius_star*radius_star);

  
  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}

}