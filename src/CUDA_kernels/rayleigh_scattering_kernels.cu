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



__global__ void rayleighScatteringCO(
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
    if (wavelength[tid] < 2.0)
    {
      const double refractive_index = (22851. + 0.456e14 / (71427.0*71427.0 - wavenumber[tid]*wavenumber[tid])) * 1e-8 + 1;

      scattering_coeff[grid_point*nb_spectral_points + tid] = 
        rayleighCrossSection(
          reference_density,
          refractive_index,
          king_correction,
          wavenumber[tid]);

      scattering_coeff[grid_point*nb_spectral_points + tid] *= number_density;
    }
  }

}



__global__ void rayleighScatteringCO2(
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

  const double reference_density = 2.546899e19; //molecules cm^-3

  if (tid < nb_spectral_points)
  {
    if (wavelength[tid] < 2.0)
    {
      const double nu = wavenumber[tid];
      const double king_correction =  1.1364 + 25.3e-12 * nu*nu;
      const double refractive_index = (5799.25 / (128908.9*128908.9 - nu*nu) 
                                    + 120.05 / (89223.8*89223.8 - nu*nu) 
                                    + 5.3334 / (75037.5*75037.5 - nu*nu) 
                                    + 4.3244 / (67837.7*67837.7 - nu*nu) 
                                    + 0.1218145e-6 / (2418.136*2418.136 - nu*nu))
                                    * 1.1427e3 + 1.0;

      scattering_coeff[grid_point*nb_spectral_points + tid] = 
        rayleighCrossSection(
          reference_density,
          refractive_index,
          king_correction,
          wavenumber[tid]);

      scattering_coeff[grid_point*nb_spectral_points + tid] *= number_density;
    }
  }

}



__global__ void rayleighScatteringCH4(
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
    if (wavelength[tid] < 2.0)
    {
      double refractive_index = (46662. + 4.02e-6 *wavenumber[tid]*wavenumber[tid]) * 1e-8 + 1;

      scattering_coeff[grid_point*nb_spectral_points + tid] = 
        rayleighCrossSection(
          reference_density,
          refractive_index,
          king_correction,
          wavenumber[tid]);

      scattering_coeff[grid_point*nb_spectral_points + tid] *= number_density;
    }
  }

}


__global__ void rayleighScatteringH2O(
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

  constexpr double king_correction = (6 + 3 * 3e-4) / (6 - 7 * 3e-4);

  //this is the number density of water at standard temperature and pressure
  //determined from the Avogradro constant and the properties of water at STP
  constexpr double reference_density = 3.34279671749673e+22; //molecules cm^-3

  //values for water at STP
  constexpr double delta = 1.0;
  constexpr double theta = 1.0;

  constexpr double lambda_uv = 0.229202;
  constexpr double lambda_ir = 5.432937;
  
  constexpr double a_coeff[8] = {0.244257733, 0.974634476e-2, -0.373234996e-2, 0.268678472e-3, 0.158920570e-2, 0.245934259e-2, 0.900704920, -0.166626219e-1};

  if (tid < nb_spectral_points)
  { 
    if (wavelength[tid] < 2.0)
    {
      const double lambda = wavelength[tid] / 0.589;

      const double a_factor = delta * (a_coeff[0] 
                            + a_coeff[1]*delta 
                            + a_coeff[2]*theta 
                            + a_coeff[3]*lambda*lambda*theta 
                            + a_coeff[4]*std::pow(lambda,-2) 
                            + a_coeff[5] / (lambda*lambda - lambda_uv*lambda_uv) 
                            + a_coeff[6] / (lambda*lambda - lambda_ir*lambda_ir) + a_coeff[7]*delta*delta);

      const double refractive_index = pow(((2 * a_factor + 1)/(1 - a_factor)), 0.5);

      scattering_coeff[grid_point*nb_spectral_points + tid] = 
        rayleighCrossSection(
          reference_density,
          refractive_index,
          king_correction,
          wavenumber[tid]);

      scattering_coeff[grid_point*nb_spectral_points + tid] *= number_density; 
    }
  }

}





__host__ void GasCO::calcRalyleighCrossSectionsGPU(
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

  rayleighScatteringCO<<<blocks,threads>>>(
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


__host__ void GasCO2::calcRalyleighCrossSectionsGPU(
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

  rayleighScatteringCO2<<<blocks,threads>>>(
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


__host__ void GasCH4::calcRalyleighCrossSectionsGPU(
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

  rayleighScatteringCH4<<<blocks,threads>>>(
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


__host__ void GasH2O::calcRalyleighCrossSectionsGPU(
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

  rayleighScatteringH2O<<<blocks,threads>>>(
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
