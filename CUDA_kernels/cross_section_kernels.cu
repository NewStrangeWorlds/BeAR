/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2019 Daniel Kitzmann
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


#include "cross_section_kernels.h"

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


//calculates the cross sections for a given p-T pair based on tabulated data
//uses a two-dimensional linar interpolation in log pressure and linear temperature based on the four closest tabulated p-T points
//note that the cross-sections are given in log10 here
__global__ void calcCrossSectionsDevice(const double* cross_sections1, const double* cross_sections2, const double* cross_sections3, const double* cross_sections4,
                                        const double temperature1, const double temperature2,
                                        const double log_pressure1, const double log_pressure2,
                                        const double temperature, const double log_pressure, const double number_density,
                                        const int nb_spectral_points, const int nb_grid_points, const int grid_point,
                                        double* __restrict__ absorption_coeff_device)
{

  //tid is the wavenumber index
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  
  if (tid < nb_spectral_points)
  {
    double c1 = cross_sections1[tid];
    double c2 = cross_sections2[tid];
    double c3 = cross_sections3[tid];
    double c4 = cross_sections4[tid];


    c1 = c1 + (c2 - c1) * (temperature - temperature1)/(temperature2 - temperature1);
    c2 = c3 + (c4 - c3) * (temperature - temperature1)/(temperature2 - temperature1);

    double sigma = c1 + (c2 - c1) * (log_pressure - log_pressure1)/(log_pressure2 - log_pressure1);
    sigma = exp10(sigma);

   
    absorption_coeff_device[grid_point*nb_spectral_points + tid] += sigma * number_density;
  }


}


//set all absorption coefficients to 0
__global__ void initCrossSectionsDevice(const int nb_points, double* absorption_coeff_device)
{

  //tid is the wavenumber index
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  
  if (tid < nb_points)
    absorption_coeff_device[tid] = 0.0;


}


//calculates the contributions due to collision induced absorption
//performs a 1D linear interpolation of the tabulated data as a function of temperature
//note that the cross-sections are given in log10 and that number_densities contains the product of both collision partners
__global__ void calcCIACoefficientsDevice(const double* cross_sections1, const double* cross_sections2,
                                          const double temperature1, const double temperature2,
                                          const double temperature, const double number_densities,
                                          const int nb_spectral_points, const int nb_grid_points, const int grid_point,
                                          double* __restrict__ absorption_coeff_device)
{

  //tid is the wavenumber index
  int tid = blockIdx.x * blockDim.x + threadIdx.x;


  if (tid < nb_spectral_points)
  {
    double coeff1 = cross_sections1[tid];
    double coeff2 = cross_sections2[tid];

    double sigma = exp10(coeff1 + (coeff2 - coeff1) * (temperature - temperature1)/(temperature2 - temperature1));


    absorption_coeff_device[grid_point*nb_spectral_points + tid] += sigma * number_densities;
  }


}



//prints the calculated cross-sections
//this is just for debug purposes
__global__ void printAC(const size_t nb_spectral_points, const size_t nb_grid_points, const size_t grid_point,
                        double* absorption_coeff_device)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  if (tid == 0)
  {


    for (int i=0; i<nb_grid_points; ++i)
    {
       //double* absorption_coeff = absorption_coeff_device[i]; //pointer to the absorption coefficient per wavenumber


       printf("%d %1.6e\n",i,absorption_coeff_device[i*nb_spectral_points + 1]);

    }


  }


}



//checks the calculated cross-sections for too large or small values
//this is just for debug purposes
__global__ void checkCrossSections(const int nb_points, double* absorption_coeff_dev)
{


  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_points; tid += blockDim.x * gridDim.x)
  {

    if (absorption_coeff_dev[tid] != absorption_coeff_dev[tid] || absorption_coeff_dev[tid] > 1e10 || absorption_coeff_dev[tid] < 0)
      printf("abs coeff are wrong: %d %1.6e\n",tid,absorption_coeff_dev[tid]);

  }


}



//sets all cross-sections back to zero
__host__ void initCrossSectionsHost(const size_t nb_points, double* absorption_coeff_device)
{


  int threads = 256;
  //int blocks = min((int (nb_points)/4+threads-1)/threads, 2048);
  int blocks = nb_points / threads;
  if (nb_points % threads) blocks++;


  initCrossSectionsDevice<<<blocks,threads>>>(nb_points, absorption_coeff_device);


  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
}



//Checks computed absorption coefficients for unrealistic values
//This is just a debug function that is normally only called to validate calculations
__host__ void checkCrossSectionsHost(const size_t nb_spectral_points, const size_t nb_grid_points,
                                     double* absorption_coeff_device)
{


  int threads = 256;
  //int blocks = min((int (nb_spectral_points)/4+threads-1)/threads, 2048);
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;


  checkCrossSections<<<blocks,threads>>>(nb_spectral_points*nb_grid_points, absorption_coeff_device);


  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
}



__host__ void calcCrossSectionsHost(double* cross_sections1, double* cross_sections2, double* cross_sections3, double* cross_sections4,
                                    const double temperature1, const double temperature2,
                                    const double pressure1, const double pressure2,
                                    const double temperature, const double pressure, const double number_density,
                                    const size_t nb_spectral_points, const size_t nb_grid_points, const size_t grid_point,
                                    double* absorption_coeff_device, double* scattering_coeff_device)

{
  double log_pressure1_gpu = log10(pressure1);
  double log_pressure2_gpu = log10(pressure2);
  double temperature1_gpu = temperature1;
  double temperature2_gpu = temperature2;

  double log_pressure = log10(pressure);
  
  //the linear interpolation in the CUDA kernel will fail if the two temperatures or pressures are equal
  //to avoid a bunch of if-statements in the CUDA kernel, we here simply offset one of the temperatures or pressures by bit
  if (temperature1_gpu == temperature2_gpu) temperature2_gpu += 1;
  if (log_pressure1_gpu == log_pressure2_gpu) log_pressure2_gpu += 0.01;


  cudaDeviceSynchronize();


  int threads = 256;
  
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;


  calcCrossSectionsDevice<<<blocks,threads>>>(cross_sections1, cross_sections2, cross_sections3, cross_sections4,
                                     temperature1_gpu, temperature2_gpu,
                                     log_pressure1_gpu, log_pressure2_gpu,
                                     temperature, log_pressure, number_density,
                                     nb_spectral_points, nb_grid_points, grid_point,
                                     absorption_coeff_device);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
}




__host__ void calcCIACoefficientsHost(double* cross_sections1, double* cross_sections2,
                                      const double temperature1, const double temperature2,
                                      const double temperature, const double number_densities,
                                      const size_t nb_spectral_points, const size_t nb_grid_points, const size_t grid_point,
                                      double* absorption_coeff_device)
{

  cudaDeviceSynchronize();


  int threads = 256;
  
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;


  double temperature1_gpu = temperature1;
  double temperature2_gpu = temperature2;

  //the linear interpolation in the CUDA kernel will fail if the two temperatures equal
  //to avoid a bunch of if-statements in the CUDA kernel, we here simply offset one of the temperatures a bit
  if (temperature2 == temperature1) temperature2_gpu += 1.0;

  calcCIACoefficientsDevice<<<blocks,threads>>>(cross_sections1, cross_sections2,
                                                temperature1_gpu, temperature2_gpu,
                                                temperature, number_densities,
                                                nb_spectral_points, nb_grid_points, grid_point,
                                                absorption_coeff_device);




  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
}






}
