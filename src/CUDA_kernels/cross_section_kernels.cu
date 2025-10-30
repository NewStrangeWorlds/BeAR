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


namespace bear{


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
    const double pressure_interpol_factor = (log_pressure - log_pressure1)/(log_pressure2 - log_pressure1);


    /*c1 = c1 + (c2 - c1) * (temperature - temperature1)/(temperature2 - temperature1);
    c2 = c3 + (c4 - c3) * (temperature - temperature1)/(temperature2 - temperature1);

    double sigma = c1 + (c2 - c1) * (log_pressure - log_pressure1)/(log_pressure2 - log_pressure1);
    sigma = exp10(sigma);*/


    c1 = c1 + (c2 - c1) * pressure_interpol_factor;
    c2 = c3 + (c4 - c3) * pressure_interpol_factor;

    double sigma = c1 + (c2 - c1) * (temperature - temperature1)/(temperature2 - temperature1);
    sigma = exp10(sigma);

   
    absorption_coeff_device[grid_point*nb_spectral_points + tid] += sigma * number_density;
  }

}



//calculates the H- continuum
__global__ void calcHmContinuumDevice(const double hm_number_density, const int nb_spectral_points, const int grid_point, double* wavelengths_d, double* __restrict__ absorption_coeff_device)
{
  const double lambda_0 = 1.6419; //photo-detachment threshold
  const double C_n[7] = {0.0, 152.519, 49.534, -118.858, 92.536, -34.194, 4.982};

  //tid is the wavenumber index
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

 
  if (tid < nb_spectral_points)
  { 
    if (wavelengths_d[tid] <= lambda_0 && wavelengths_d[tid] >= 0.125)
    {
      const double lambda = wavelengths_d[tid];
      double f = 0;

      const double x = 1.0/lambda - 1.0/lambda_0;

      for (unsigned int i=1; i<7; ++i)
        f += C_n[i] * pow(x, (i-1)/2.0);

      const double sigma = 1e-18 * wavelengths_d[tid] * wavelengths_d[tid] * wavelengths_d[tid] * std::pow(x, 1.5) * f; 
      
      absorption_coeff_device[grid_point*nb_spectral_points + tid] += sigma * hm_number_density;
    }

  }


}



//calculates the free-free H- continuum
__global__ void calcHmffContinuumDevice(const double h_number_density, const double e_pressure, const double temperature,
                                        const int nb_spectral_points, const int grid_point, double* wavelengths_d, double* __restrict__ absorption_coeff_device)
{
  //for wavelengths larger than 0.3645 micron
  const double A_n1[7] = {0.0, 0.0, 2483.3460, -3449.8890, 2200.0400, -696.2710, 88.2830};
  const double B_n1[7] = {0.0, 0.0, 285.8270, -1158.3820, 2427.7190, -1841.4000, 444.5170};
  const double C_n1[7] = {0.0, 0.0, -2054.2910, 8746.5230, -13651.1050, 8624.9700, -1863.8650};
  const double D_n1[7] = {0.0, 0.0, 2827.7760, -11485.6320, 16755.5240, -10051.5300, 2095.2880};
  const double E_n1[7] = {0.0, 0.0, -1341.5370, 5303.6090, -7510.4940, 4400.0670, -901.7880};
  const double F_n1[7] = {0.0, 0.0, 208.9520, -812.9390, 1132.7380, -655.0200, 132.9850};

  //for wavelengths between 0.1823 micron and 0.3645 micron
  const double A_n2[7] = {0.0, 518.1021, 473.2636, -482.2089, 115.5291, 0.0, 0.0};
  const double B_n2[7] = {0.0, -734.8666, 1443.4137, -737.1616, 169.6374, 0.0, 0.0};
  const double C_n2[7] = {0.0, 1021.1775, -1977.3395, 1096.8827, -245.6490, 0.0, 0.0};
  const double D_n2[7] = {0.0, -479.0721, 922.3575, -521.1341, 114.2430, 0.0, 0.0};
  const double E_n2[7] = {0.0, 93.1373, -178.9275, 101.7963, -21.9972, 0.0, 0.0};
  const double F_n2[7] = {0.0, -6.4285, 12.3600, -7.0571, 1.5097, 0.0, 0.0};


  //tid is the wavenumber index
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  if (tid < nb_spectral_points)
  {
    const double lambda = wavelengths_d[tid];

    double ff_sigma = 0;

    if (lambda > 0.3645)
    {
      for (unsigned int i=1; i<7; ++i)
        ff_sigma += std::pow(5040.0/temperature, (i+1)/2.0) * 
             (lambda * lambda * A_n1[i] + B_n1[i] + C_n1[i]/lambda + D_n1[i]/lambda/lambda 
               + E_n1[i]/lambda/lambda/lambda + F_n1[i]/lambda/lambda/lambda/lambda);
    }
    else if (lambda >= 0.1823)
    {
      for (unsigned int i=1; i<7; ++i)
        ff_sigma += std::pow(5040.0/temperature, (i+1)/2.0) * 
              (lambda * lambda * A_n2[i] + B_n2[i] + C_n2[i]/lambda + D_n2[i]/lambda/lambda 
              + E_n2[i]/lambda/lambda/lambda + F_n2[i]/lambda/lambda/lambda/lambda);
    }
    
    ff_sigma *= 1e-29;
    if (ff_sigma < 0) ff_sigma = 0;
   
    absorption_coeff_device[grid_point*nb_spectral_points + tid] += ff_sigma * h_number_density * e_pressure;
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
                                    const double log_pressure1, const double log_pressure2,
                                    const double temperature, const double log_pressure, const double number_density,
                                    const size_t nb_spectral_points, const size_t nb_grid_points, const size_t grid_point,
                                    double* absorption_coeff_device, double* scattering_coeff_device)

{
  double log_pressure1_gpu = log_pressure1;
  double log_pressure2_gpu = log_pressure2;
  double temperature1_gpu = temperature1;
  double temperature2_gpu = temperature2;

  
  //the linear interpolation in the CUDA kernel will fail if the two temperatures or pressures are equal
  //to avoid a bunch of if-statements in the CUDA kernel, we here simply offset one of the temperatures or pressures by a bit
  if (temperature1_gpu == temperature2_gpu) temperature2_gpu += 1;
  if (log_pressure1_gpu == log_pressure2_gpu) log_pressure2_gpu += 0.001;


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



__host__ void calcHmContinuumHost(const double hm_number_density,
                                  const double h_number_density,
                                  const double e_pressure,
                                  const double temperature,
                                  const int nb_spectral_points, const int grid_point,
                                  double* wavelengths_device,
                                  double* absorption_coeff_device)
{

int threads = 256;
//int blocks = min((int (nb_spectral_points)/4+threads-1)/threads, 2048);
int blocks = nb_spectral_points / threads;
if (nb_spectral_points % threads) blocks++;


calcHmContinuumDevice<<<blocks,threads>>>(hm_number_density, nb_spectral_points, grid_point, wavelengths_device, absorption_coeff_device);


cudaDeviceSynchronize();

calcHmffContinuumDevice<<<blocks,threads>>>(h_number_density, e_pressure, temperature, nb_spectral_points, grid_point, wavelengths_device, absorption_coeff_device);


cudaDeviceSynchronize();
gpuErrchk( cudaPeekAtLastError() );
}







}
