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
#include <cmath>
#include <stdio.h>

#include "../forward_model/secondary_eclipse_bb/secondary_eclipse_bb.h"

#include "error_check.h"
#include "planck_function.h"
#include "../additional/physical_const.h"


namespace bear{


__global__ void planetBlackBodyFlux(
  double* wavenumbers,
  int nb_wavenumbers,
  const double planet_temperature,
  double* flux)
{
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_wavenumbers; tid += blockDim.x * gridDim.x)
  {
    flux[tid] = planckFunction(
      planet_temperature,
      wavenumbers[tid]) * constants::pi * 1e-3;
  }
}



__global__ void secondaryEclipseBBDevice(
  double* secondary_eclipse,
  double* planet_spectrum,
  const double* stellar_spectrum,
  const int nb_points,
  const double radius_ratio, 
  const double* albedo_contribution)
{
  
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < nb_points; i += blockDim.x * gridDim.x)
  {
    secondary_eclipse[i] = planet_spectrum[i]/stellar_spectrum[i] * radius_ratio*radius_ratio * 1e6; 

    if (albedo_contribution != nullptr)
      secondary_eclipse[i] += albedo_contribution[i]*1e6;
  }
}


__host__ void SecondaryEclipseBlackBodyModel::calcPlanetSpectrumGPU(
  const double planet_temperature,
  double* spectrum_dev)
{
  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  planetBlackBodyFlux<<<blocks,threads>>>(
    spectral_grid->wavenumber_list_gpu,
    nb_spectral_points,
    planet_temperature,
    spectrum_dev);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}


__host__ void SecondaryEclipseBlackBodyModel::calcSecondaryEclipseGPU(
  double* secondary_eclipse,
  double* planet_spectrum,
  const double* stellar_spectrum,
  const int nb_points,
  const double radius_ratio,
  const double* albedo_contribution)
{
  int threads = 256;

  int blocks = nb_points / threads;
  if (nb_points % threads) blocks++;


  secondaryEclipseBBDevice<<<blocks,threads>>>(
    secondary_eclipse,
    planet_spectrum,
    stellar_spectrum,
    nb_points,
    radius_ratio,
    albedo_contribution);


  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}


}