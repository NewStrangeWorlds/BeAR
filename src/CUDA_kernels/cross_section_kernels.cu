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
#include "../transport_coeff/species_definition.h"
#include "../transport_coeff/opacity_species.h"


namespace bear{


// H- bound-free coefficients (constant memory)
__constant__ double d_C_n[7] = {0.0, 152.519, 49.534, -118.858, 92.536, -34.194, 4.982};

// H- free-free coefficients for lambda > 0.3645 micron (constant memory)
__constant__ double d_A_n1[7] = {0.0, 0.0, 2483.3460, -3449.8890, 2200.0400, -696.2710, 88.2830};
__constant__ double d_B_n1[7] = {0.0, 0.0, 285.8270, -1158.3820, 2427.7190, -1841.4000, 444.5170};
__constant__ double d_C_n1[7] = {0.0, 0.0, -2054.2910, 8746.5230, -13651.1050, 8624.9700, -1863.8650};
__constant__ double d_D_n1[7] = {0.0, 0.0, 2827.7760, -11485.6320, 16755.5240, -10051.5300, 2095.2880};
__constant__ double d_E_n1[7] = {0.0, 0.0, -1341.5370, 5303.6090, -7510.4940, 4400.0670, -901.7880};
__constant__ double d_F_n1[7] = {0.0, 0.0, 208.9520, -812.9390, 1132.7380, -655.0200, 132.9850};

// H- free-free coefficients for 0.1823 <= lambda <= 0.3645 micron (constant memory)
__constant__ double d_A_n2[7] = {0.0, 518.1021, 473.2636, -482.2089, 115.5291, 0.0, 0.0};
__constant__ double d_B_n2[7] = {0.0, -734.8666, 1443.4137, -737.1616, 169.6374, 0.0, 0.0};
__constant__ double d_C_n2[7] = {0.0, 1021.1775, -1977.3395, 1096.8827, -245.6490, 0.0, 0.0};
__constant__ double d_D_n2[7] = {0.0, -479.0721, 922.3575, -521.1341, 114.2430, 0.0, 0.0};
__constant__ double d_E_n2[7] = {0.0, 93.1373, -178.9275, 101.7963, -21.9972, 0.0, 0.0};
__constant__ double d_F_n2[7] = {0.0, -6.4285, 12.3600, -7.0571, 1.5097, 0.0, 0.0};


//calculates the cross sections for a given p-T pair based on tabulated data
//uses a two-dimensional linar interpolation in log pressure and linear temperature based on the four closest tabulated p-T points
//note that the cross-sections are given in log10 here
__global__ void calcCrossSectionsDevice(
  const float* __restrict__ cross_sections1, 
  const float* __restrict__ cross_sections2, 
  const float* __restrict__ cross_sections3, 
  const float* __restrict__ cross_sections4,
  const float temperature_interpol_factor,
  const float pressure_interpol_factor,
  const double number_density,
  const int nb_spectral_points, 
  const int grid_point,
  float* __restrict__ absorption_coeff_device)
{
  //tid is the wavenumber index
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  if (tid < nb_spectral_points)
  {
    float c1 = cross_sections1[tid];
    float c2 = cross_sections2[tid];
    const float c3 = cross_sections3[tid];
    const float c4 = cross_sections4[tid];
    
    c1 = c1 + (c2 - c1) * pressure_interpol_factor;
    c2 = c3 + (c4 - c3) * pressure_interpol_factor;

    double sigma = c1 + (c2 - c1) * temperature_interpol_factor;
    
    sigma = exp10(sigma) * number_density;
   
    absorption_coeff_device[grid_point*nb_spectral_points + tid] += sigma;
  }
}


//calculates the H- continuum
__global__ void calcHmbfContinuumDevice(
  const double hm_number_density,
  const int nb_spectral_points,
  const int grid_point,
  float* __restrict__ cross_section_dev,
  float* __restrict__ absorption_coeff_device)
{
  //tid is the wavenumber index
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  if (tid < nb_spectral_points)
  {
    absorption_coeff_device[grid_point*nb_spectral_points + tid] += 
      cross_section_dev[tid] * hm_number_density;
  }
}


//tabulate the H- bound-free cross-section based on the formula from John (1988) 
//and the coefficients from Bell & Berrington (1987)
__global__ void calcHmbfCrossSectionDevice(
  const int nb_spectral_points,
  double* wavelengths_dev,
  float* __restrict__ cross_section_dev)
{
  const double lambda_0 = 1.6419; //photo-detachment threshold

  //tid is the wavenumber index
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  if (tid < nb_spectral_points)
  {
    const double lambda = wavelengths_dev[tid];

    if (lambda <= lambda_0 && lambda >= 0.125)
    {
      const double x = 1.0/lambda - 1.0/lambda_0;
      const double sqrt_x = sqrt(x);

      //powers of sqrt(x): x^0, x^0.5, x^1, x^1.5, x^2, x^2.5
      //computed incrementally instead of calling pow() 6 times
      //f = C_n[1]*x^0 + C_n[2]*x^0.5 + C_n[3]*x^1 + C_n[4]*x^1.5 + C_n[5]*x^2 + C_n[6]*x^2.5
      double x_power = 1.0;  // x^0
      double f = d_C_n[1];   // i=1: C_n[1] * x^0

      x_power *= sqrt_x;     // x^0.5
      f += d_C_n[2] * x_power;

      x_power *= sqrt_x;     // x^1
      f += d_C_n[3] * x_power;

      x_power *= sqrt_x;     // x^1.5
      f += d_C_n[4] * x_power;

      x_power *= sqrt_x;     // x^2
      f += d_C_n[5] * x_power;

      x_power *= sqrt_x;     // x^2.5
      f += d_C_n[6] * x_power;

      const double lambda3 = lambda * lambda * lambda;
      const double x15 = x * sqrt_x;  // x^1.5
      const double sigma = 1e-18 * lambda3 * x15 * f;

      cross_section_dev[tid] = sigma;
    }
  }
}



//calculates the free-free H- continuum
__global__ void calcHmffContinuumDevice(
  const double h_number_density,
  const double e_pressure,
  const double temperature,
  const int nb_spectral_points,
  const int grid_point,
  double* wavelengths_d,
  float* __restrict__ absorption_coeff_device)
{
  //tid is the wavenumber index
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  if (tid < nb_spectral_points)
  {
    const double lambda = wavelengths_d[tid];

    double ff_sigma = 0;

    if (lambda >= 0.1823)
    {
      //precompute inverse lambda powers once instead of repeated divisions
      const double inv_l = 1.0 / lambda;
      const double l2 = lambda * lambda;
      const double inv_l2 = inv_l * inv_l;
      const double inv_l3 = inv_l2 * inv_l;
      const double inv_l4 = inv_l2 * inv_l2;

      //precompute theta powers incrementally: theta^1, theta^1.5, theta^2, ..., theta^3.5
      //instead of calling pow() 6 times
      const double theta = 5040.0 / temperature;
      const double sqrt_theta = sqrt(theta);
      double theta_power = theta;  // theta^1 for i=1: (i+1)/2 = 1

      //select correct coefficient set based on wavelength
      const bool use_set1 = (lambda > 0.3645);

      #pragma unroll
      for (unsigned int i=1; i<7; ++i)
      {
        const double A_i = use_set1 ? d_A_n1[i] : d_A_n2[i];
        const double B_i = use_set1 ? d_B_n1[i] : d_B_n2[i];
        const double C_i = use_set1 ? d_C_n1[i] : d_C_n2[i];
        const double D_i = use_set1 ? d_D_n1[i] : d_D_n2[i];
        const double E_i = use_set1 ? d_E_n1[i] : d_E_n2[i];
        const double F_i = use_set1 ? d_F_n1[i] : d_F_n2[i];

        const double lambda_poly = l2 * A_i + B_i + C_i * inv_l + D_i * inv_l2
                                   + E_i * inv_l3 + F_i * inv_l4;

        ff_sigma += theta_power * lambda_poly;
        theta_power *= sqrt_theta;  // next power: theta^1.5, theta^2, ...
      }

      ff_sigma *= 1e-29;
      if (ff_sigma < 0) ff_sigma = 0;

      absorption_coeff_device[grid_point*nb_spectral_points + tid] += ff_sigma * h_number_density * e_pressure;
    }
  }
}



//set all absorption coefficients to 0
__global__ void initCrossSectionsDevice(
  const int nb_points, 
  float* absorption_coeff_device)
{
  //tid is the wavenumber index
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  if (tid < nb_points)
    absorption_coeff_device[tid] = 0.0;
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
__global__ 
void checkCrossSections(
  const int nb_points, 
  float* absorption_coeff_dev)
{
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_points; tid += blockDim.x * gridDim.x)
  {
    if (absorption_coeff_dev[tid] != absorption_coeff_dev[tid] || absorption_coeff_dev[tid] > 1e10 || absorption_coeff_dev[tid] < 0)
      printf("abs coeff are wrong: %d %1.6e\n",tid,absorption_coeff_dev[tid]);
  }

}



//sets all cross-sections back to zero
__host__ void initCrossSectionsHost(
  const size_t nb_points, 
  float* absorption_coeff_device)
{
  cudaMemset(absorption_coeff_device, 0, nb_points * sizeof(float));

  gpuErrchk( cudaPeekAtLastError() );
}



//Checks computed absorption coefficients for unrealistic values
//This is just a debug function that is normally only called to validate calculations
__host__ void checkCrossSectionsHost(
  const size_t nb_spectral_points, 
  const size_t nb_grid_points,
  float* absorption_coeff_device)
{
  int threads = 256;
  //int blocks = min((int (nb_spectral_points)/4+threads-1)/threads, 2048);
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  checkCrossSections<<<blocks,threads>>>(nb_spectral_points*nb_grid_points, absorption_coeff_device);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
}



__host__ void OpacitySpecies::calcAbsorptionCoefficientsGPU(
  float* cross_sections1, 
  float* cross_sections2, 
  float* cross_sections3, 
  float* cross_sections4,
  const double temperature1, 
  const double temperature2,
  const double log_pressure1, 
  const double log_pressure2,
  const double temperature, 
  const double log_pressure, 
  const double number_density,
  const size_t nb_spectral_points, 
  const size_t nb_grid_points, 
  const size_t grid_point,
  float* absorption_coeff_device)
{
  double log_pressure1_gpu = log_pressure1;
  double log_pressure2_gpu = log_pressure2;
  double temperature1_gpu = temperature1;
  double temperature2_gpu = temperature2;

  //the linear interpolation in the CUDA kernel will fail if the two temperatures 
  //or pressures are equal to avoid a bunch of if-statements in the CUDA kernel, 
  //we here simply offset one of the temperatures or pressures by a bit
  if (temperature1_gpu == temperature2_gpu) temperature2_gpu += 1;
  if (log_pressure1_gpu == log_pressure2_gpu) log_pressure2_gpu += 0.001;

  int threads = 256;
  
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  const double pressure_interpol_factor = 
    (log_pressure - log_pressure1_gpu)/(log_pressure2_gpu - log_pressure1_gpu);
  const double temperature_interpol_factor = 
    (temperature - temperature1_gpu)/(temperature2_gpu - temperature1_gpu);


  calcCrossSectionsDevice<<<blocks,threads>>>(
    cross_sections1, 
    cross_sections2, 
    cross_sections3, 
    cross_sections4,
    temperature_interpol_factor,
    pressure_interpol_factor,
    number_density,
    nb_spectral_points,
    grid_point,
    absorption_coeff_device);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
}




__host__ void GasHm::calcContinuumGPU(
  const double hm_number_density,
  const double h_number_density,
  const double e_pressure,
  const double temperature,
  const int nb_spectral_points, 
  const int grid_point,
  double* wavelengths_device,
  float* absorption_coeff_device)
{
  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;
  
  //pre-tabulate the bound-free cross-section if not already done
  if (bound_free_sigma_dev == nullptr)
  {
    allocateOnDevice(bound_free_sigma_dev, nb_spectral_points * sizeof(float));

    calcHmbfCrossSectionDevice<<<blocks,threads>>>(
      nb_spectral_points,
      wavelengths_device,
      bound_free_sigma_dev);
  }

  calcHmbfContinuumDevice<<<blocks,threads>>>(
    hm_number_density,
    nb_spectral_points,
    grid_point,
    bound_free_sigma_dev,
    absorption_coeff_device);

  calcHmffContinuumDevice<<<blocks,threads>>>(
    h_number_density, 
    e_pressure, 
    temperature,
    nb_spectral_points, 
    grid_point,
    wavelengths_device, 
    absorption_coeff_device);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
}







}
