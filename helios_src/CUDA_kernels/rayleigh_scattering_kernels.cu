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


namespace helios{


__forceinline__ __device__ double rayleighCrossSection(
  const double reference_density, 
  const double refractive_index, 
  const double king_correction,
  const double wavenumber)
{
  return 24. * pow(constants::pi, 3) * pow(wavenumber, 4) / (reference_density*reference_density)
         * pow((refractive_index*refractive_index - 1) / (refractive_index*refractive_index + 2),2) * king_correction;
}



__global__ void rayleighScatteringH2(
  const double number_density,
  const int nb_spectral_points,
  const int nb_grid_points, 
  const int grid_point,
  double* wavelength,
  double* wavenumber,
  double* scattering_coeff)
{
  //tid is the wavenumber index
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  
  const double king_correction = 1.0;
  const double reference_density = 2.651629e19; //molecules cm^-3

  if (tid < nb_spectral_points)
  { 
    const double refractive_index = (13.58e-5 * (1. + 7.52e-3/std::pow(wavelength[tid],2))) + 1.;
    
    scattering_coeff[grid_point*nb_spectral_points + tid] = 
      rayleighCrossSection(
        reference_density,
        refractive_index,
        king_correction,
        wavenumber[tid]);
    
    scattering_coeff[grid_point*nb_spectral_points + tid] *= number_density; 
  }

}



__global__ void rayleighScatteringHe(
  const double number_density,
  const int nb_spectral_points,
  const int nb_grid_points, 
  const int grid_point,
  double* wavelength,
  double* wavenumber,
  double* scattering_coeff)
{
  //tid is the wavenumber index
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  const double king_correction = 1.0;
  const double reference_density = 2.546899e19; //molecules cm^-3

  if (tid < nb_spectral_points)
  { 
    const double refractive_index = (2283. + 1.8102e13 / (1.5342e10 - wavenumber[tid]*wavenumber[tid])) * 1e-8 + 1;

    scattering_coeff[grid_point*nb_spectral_points + tid] = 
      rayleighCrossSection(
        reference_density,
        refractive_index,
        king_correction,
        wavenumber[tid]);

    scattering_coeff[grid_point*nb_spectral_points + tid] *= number_density; 
  }

}



__host__ void GasH2::calcRalyleighCrossSectionsGPU(
  const double number_density,
  const size_t nb_grid_points, 
  const size_t grid_point,
  double* scattering_coeff)
{
  cudaDeviceSynchronize();

  const auto nb_spectral_points = spectral_grid->nbSpectralPoints();
  int threads = 256;
  
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  rayleighScatteringH2<<<blocks,threads>>>(
    number_density,
    nb_spectral_points,
    nb_grid_points,
    grid_point,
    spectral_grid->wavelength_list_gpu,
    spectral_grid->wavenumber_list_gpu,
    scattering_coeff);


  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
}



__host__ void GasHe::calcRalyleighCrossSectionsGPU(
  const double number_density,
  const size_t nb_grid_points, 
  const size_t grid_point,
  double* scattering_coeff)
{
  cudaDeviceSynchronize();

  const auto nb_spectral_points = spectral_grid->nbSpectralPoints();
  int threads = 256;
  
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  rayleighScatteringHe<<<blocks,threads>>>(
    number_density,
    nb_spectral_points,
    nb_grid_points,
    grid_point,
    spectral_grid->wavelength_list_gpu,
    spectral_grid->wavenumber_list_gpu,
    scattering_coeff);


  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
}




}
