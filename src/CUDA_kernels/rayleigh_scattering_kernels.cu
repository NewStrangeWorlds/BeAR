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


#include "../transport_coeff/species_definition.h"

#include <iostream>
#include <vector>
#include "math.h"
#include <stdio.h>
#include <new>

#include <fstream>
#include <string>
#include <iomanip>

#include "error_check.h"
#include "../additional/physical_const.h"


namespace bear{


__global__ void rayleighScatteringDevice(
  const double number_density,
  const int nb_spectral_points,
  const int grid_point,
  double* rayleigh_cross_sections_dev,
  double* scattering_coeff)
{
  //tid is the wavenumber index
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  if (tid < nb_spectral_points)
  { 
    scattering_coeff[grid_point*nb_spectral_points + tid] += rayleigh_cross_sections_dev[tid] * number_density; 
  }

}



__host__ void  OpacitySpecies::calcRayleighCrossSectionsGPU(
  const double number_density,
  const size_t nb_grid_points, 
  const size_t grid_point,
  double* scattering_coeff_dev)
{
  if (rayleigh_available == false)
    return;
  
  cudaDeviceSynchronize();

  const auto nb_spectral_points = spectral_grid->nbSpectralPoints();
  int threads = 256;
  
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  rayleighScatteringDevice<<<blocks,threads>>>(
    number_density,
    nb_spectral_points,
    grid_point,
    rayleigh_cross_sections_dev,
    scattering_coeff_dev);


  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
}



}
