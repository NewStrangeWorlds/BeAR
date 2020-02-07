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


#include "convolution_kernels.h"


#include <iostream>
#include <vector>
#include "math.h"
#include <stdio.h>
#include <new>

#include "../additional/physical_const.h"


#include "error_check.h"
#include "reduce_kernels.h"


namespace helios{


__forceinline__ __device__ double normalDistribution(const double sigma, const double x)
{
  
  const double c = sqrt(1.0/(2.0 * constants::pi));
  
  return c / sigma * exp(- x*x / (2.0 * sigma*sigma));
  
}



//each block convolves one specific wavelength
__global__ void convolveSpectrumDevice(double* spectrum, 
                                       double* band_wavelengths, double* band_sigma,
                                       int* band_indices, int* start_index, int* end_index,
                                       double* convolved_spectrum)
{
  //the current wavelength index
  const int i = blockIdx.x;

  //the current wavelength
  const double mu = band_wavelengths[i];
  
  const double sigma = band_sigma[i];

  //the length of the spectrum we need to convolve
  //we only integrate a part of the full spectrum (up to a certain number of sigmas from the central wavelength)
  const int start = start_index[i];
  const int end = end_index[i];
  const int sub_spectrum_size = end - start + 1; 

  
  //the pointer to the data vector shared by all threads
  __shared__ double* data;

  //first thread in the block does the allocation
  if (threadIdx.x == 0)
  {
    data = (double*) malloc(sub_spectrum_size * sizeof(double));
  

    //this probably needs a more sophisticated backup procedure
    if (data == nullptr)
    {
      printf("Not enough memory on GPU! %d %d %d %lu\n", start, end, sub_spectrum_size, sub_spectrum_size * sizeof(double));
      return;
    }
    
  } 
    

  __syncthreads();


  //we create the data vector for the convolution
  for (int j = threadIdx.x; j < sub_spectrum_size; j += blockDim.x)
  {
    //the wavelength index of the full spectrum
    const int index = band_indices[start + j];

    //distance from the central wavelength
    const double distance = abs(mu - band_wavelengths[start + j]);

    //the data for the convolution 
    data[j] = spectrum[index] * normalDistribution(sigma, distance);
  }


  __syncthreads();

  
  //now we integrate with a trapezoidal rule
  double quad_sum = 0;

  for (int j = threadIdx.x; j < sub_spectrum_size-1; j += blockDim.x) 
    quad_sum += (data[j] + data[j+1]) * (band_wavelengths[start + j + 1] - band_wavelengths[start + j]);
  

  __syncthreads();


  quad_sum = blockReduceSum(quad_sum);


  if (threadIdx.x == 0)
    convolved_spectrum[band_indices[i]] = abs(quad_sum * 0.5);


  //first thread frees the memory
  if (threadIdx.x == 0)
    free(data);
}






__host__ void convolveSpectrumGPU(double* spectrum, 
                                  double* band_wavelengths, double* band_sigma,
                                  int* band_indices, int* start_index, int* end_index,
                                  const int nb_points,
                                  double* convolved_spectrum)
{
  int threads = 128;
  int blocks = nb_points;


  convolveSpectrumDevice<<<blocks,threads>>>(spectrum, band_wavelengths, band_sigma,
                                             band_indices, start_index, end_index,
                                             convolved_spectrum);


  cudaDeviceSynchronize(); 
  gpuErrchk( cudaPeekAtLastError() ); 
}






}