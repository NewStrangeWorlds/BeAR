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
#include <vector>
#include "math.h"
#include <stdio.h>
#include <new>

#include "../spectral_grid/spectral_band.h"
#include "../additional/physical_const.h"


#include "error_check.h"
#include "reduce_kernels.h"
#include "../spectral_grid/spectral_grid.h"


namespace bear{


__forceinline__ __device__ double normalDistribution(const double sigma, const double x)
{
  const double c = sqrt(1.0/(2.0 * constants::pi));
  
  return c / sigma * exp(- x*x / (2.0 * sigma*sigma));
}



//each block convolves one specific wavelength
__global__ void convolveSpectrumDevice(
  double* spectrum, 
  const int index_start,
  double* wavelengths, 
  double* band_sigma,
  int* start_index, 
  int* end_index,
  double* convolved_spectrum)
{ 
  //the current wavelength index
  const int i = blockIdx.x;

  //the current wavelength
  const double mu = wavelengths[i + index_start];
  
  const double sigma = band_sigma[i + index_start];

  //the length of the spectrum we need to convolve
  //we only integrate a part of the full spectrum (up to a certain number of sigmas from the central wavelength)
  const int start = start_index[i];
  const int end = end_index[i];
  const int sub_spectrum_size = end - start + 1;
  
  if (start == end || sigma == 0)
  {
    if (threadIdx.x == 0)
      convolved_spectrum[i+index_start] = spectrum[i+index_start];

    __syncthreads();

    return;
  }


  //the pointer to the data vector shared by all threads
  __shared__ double* data;

  //first thread in the block does the allocation
  if (threadIdx.x == 0)
  {
    data = (double*) malloc(sub_spectrum_size * sizeof(double));

    //this probably needs a more sophisticated backup procedure
    if (data == nullptr)
    {
      printf("Not enough memory on GPU! %d %d %d %lu\n", 
        start, 
        end, 
        sub_spectrum_size, 
        sub_spectrum_size * sizeof(double));
      
      return;
    }
    
  } 


  __syncthreads();


  //we create the data vector for the convolution
  for (int j = threadIdx.x; j < sub_spectrum_size; j += blockDim.x)
  {
    //distance from the central wavelength
    const double distance = abs(mu - wavelengths[j + start]);

    //the data for the convolution 
    data[j] = spectrum[j+start] * normalDistribution(sigma, distance);
  }

  
  __syncthreads();

  //now we integrate with a trapezoidal rule
  double quad_sum = 0;

  for (int j = threadIdx.x; j < sub_spectrum_size-1; j += blockDim.x) 
    quad_sum += (data[j] + data[j+1]) * (wavelengths[j + 1 + start] - wavelengths[j + start]);
  

  __syncthreads();


  quad_sum = blockReduceSum(quad_sum);


  if (threadIdx.x == 0)
    convolved_spectrum[i + index_start] = abs(quad_sum * 0.5);


  //first thread frees the memory
  if (threadIdx.x == 0)
    free(data);
}



__host__ void SpectralBands::convolveSpectrumGPU(double* spectrum, double* spectrum_processed_dev)
{
  const size_t nb_high_res_points = obs_index_range.second - obs_index_range.first + 1;

  int threads = 256; //128;
  int blocks = nb_high_res_points;

  convolveSpectrumDevice<<<blocks,threads>>>(
    spectrum, 
    obs_index_range.first,
    spectral_grid->wavelength_list_gpu, 
    instrument_profile_sigma_dev,
    convolution_start_dev, 
    convolution_end_dev,
    spectrum_processed_dev);


  cudaDeviceSynchronize(); 
  gpuErrchk( cudaPeekAtLastError() ); 
}



}
