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


#include "../spectral_grid/spectral_band.h"
#include "../spectral_grid/spectral_grid.h"


#include <iostream>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <new>

#include "data_management_kernels.h"
#include "error_check.h"
#include "reduce_kernels.h"



namespace bear{

//every block reduces one band
//spectrum_high_res is the pointer to the array on the GPU
//the integration is done by a simple piece-wise trapezoidal rule
//note that the units of the high-res spectrum are in W m-2 cm, while the mean band values are in W m-2 mu-1
__global__ 
void bandIntegrationDeviceOld(
  const double* __restrict__ spectrum_high_res, 
  const int* __restrict__ band_start, 
  const int*  __restrict__ band_end,
  const double* __restrict__ wavenumbers, 
  const double* __restrict__ wavelengths,
  double* __restrict__ spectrum_bands,
  const bool is_flux,
  const bool use_filter_transmission)
{
  float band_sum = 0;
  
  //indices to navigate through the high-res spectrum
  const int start_index = band_start[blockIdx.x];
  const int end_index = band_end[blockIdx.x];
  const int band_size = end_index - start_index + 1;

  for (int j = threadIdx.x; j < band_size-1; j += blockDim.x)
  {
    const int index1 = j + start_index + 1;
    const int index2 = j + start_index;

    float delta = 0;

    if (is_flux) 
      delta = (wavenumbers[index1] - wavenumbers[index2]);
    else
      delta = (wavelengths[index2] - wavelengths[index1]);

    const float sum = (spectrum_high_res[index1] + spectrum_high_res[index2]) * delta;

    band_sum += sum;
  }

  __syncthreads();

  band_sum = blockReduceSum(band_sum);

  if (threadIdx.x == 0)
  {
    if (is_flux)
      spectrum_bands[blockIdx.x] = band_sum * 0.5 / (wavelengths[start_index] - wavelengths[end_index]);
    else 
    {
      if (use_filter_transmission)
        spectrum_bands[blockIdx.x] = band_sum * 0.5;
      else
        spectrum_bands[blockIdx.x] = band_sum * 0.5 / (wavelengths[start_index] - wavelengths[end_index]);

    }
      
  }

}


//every block reduces one band
//spectrum_high_res is the pointer to the array on the GPU
//the integration is done by a simple piece-wise trapezoidal rule
//note that the units of the high-res spectrum are in W m-2 cm, while the mean band values are in W m-2 mu-1
__global__ 
void bandIntegrationDeviceFlux(
  const double* __restrict__ spectrum_high_res, 
  const int* __restrict__ band_start, 
  const int*  __restrict__ band_end,
  const double* __restrict__ wavenumbers,
  const double* __restrict__ wavelengths,
  double* __restrict__ spectrum_bands)
{
  float band_sum = 0;
  
  //indices to navigate through the high-res spectrum
  const int start_index = band_start[blockIdx.x];
  const int end_index = band_end[blockIdx.x];
  const int band_size = end_index - start_index + 1;

  for (int j = threadIdx.x; j < band_size-1; j += blockDim.x)
  {
    const int index1 = j + start_index + 1;
    const int index2 = j + start_index;

    const float w1 = wavenumbers[index1];
    const float w2 = wavenumbers[index2];

    const float spectrum1 = spectrum_high_res[index1];
    const float spectrum2 = spectrum_high_res[index2];
    
    const float delta = (w1 - w2);
    
    const float sum = (spectrum1 + spectrum2) * delta;
    band_sum += sum;
  }

  band_sum = blockReduceSum(band_sum);

  if (threadIdx.x == 0)
  { 
    const float w1 = wavelengths[start_index];
    const float w2 = wavelengths[end_index];
    
    spectrum_bands[blockIdx.x] = band_sum * 0.5 / (w1 - w2);
  }

}



//every block reduces one band
//spectrum_high_res is the pointer to the array on the GPU
//the integration is done by a simple piece-wise trapezoidal rule
//note that the units of the high-res spectrum are in W m-2 cm, while the mean band values are in W m-2 mu-1
__global__ 
void bandIntegrationDevice(
  const double* __restrict__ spectrum_high_res, 
  const int* __restrict__ band_start, 
  const int*  __restrict__ band_end,
  const double* __restrict__ wavelengths,
  double* __restrict__ spectrum_bands,
  const bool use_filter_transmission)
{
  float band_sum = 0;
  
  //indices to navigate through the high-res spectrum
  const int start_index = band_start[blockIdx.x];
  const int end_index = band_end[blockIdx.x];
  const int band_size = end_index - start_index + 1;

  for (int j = threadIdx.x; j < band_size-1; j += blockDim.x)
  {
    const int index1 = j + start_index + 1;
    const int index2 = j + start_index;

    const float w1 = wavelengths[index1];
    const float w2 = wavelengths[index2];

    const float spectrum1 = spectrum_high_res[index1];
    const float spectrum2 = spectrum_high_res[index2];

    const float delta = (w2 - w1);
    
    const float sum = (spectrum1 + spectrum2) * delta;
    band_sum += sum;
  }

  band_sum = blockReduceSum(band_sum);

  if (threadIdx.x == 0)
  {
    if (use_filter_transmission)
      spectrum_bands[blockIdx.x] = band_sum * 0.5;
    else
    {
      const float w1 = wavelengths[start_index];
      const float w2 = wavelengths[end_index];

      spectrum_bands[blockIdx.x] = band_sum * 0.5 / (w1 - w2);
    }
  }

}




__host__ 
void SpectralBands::bandIntegrateSpectrumGPU(
  double* spectrum, 
  double* spectrum_bands, 
  const bool is_flux,
  const bool use_filter_transmission)
{
  int threads = 128;
  int blocks = nb_bands;

  // bandIntegrationDeviceOld<<<blocks,threads>>>(
  //   spectrum, 
  //   band_start_dev, 
  //   band_end_dev,
  //   spectral_grid->wavenumber_list_gpu,
  //   spectral_grid->wavelength_list_gpu,
  //   spectrum_bands,
  //   is_flux,
  //   use_filter_transmission);
  
  if (is_flux)
  {
    bandIntegrationDeviceFlux<<<blocks,threads>>>(
      spectrum, 
      band_start_dev, 
      band_end_dev,
      spectral_grid->wavenumber_list_gpu,
      spectral_grid->wavelength_list_gpu,
      spectrum_bands);
  }
  else
  {
    bandIntegrationDevice<<<blocks,threads>>>(
      spectrum, 
      band_start_dev, 
      band_end_dev,
      spectral_grid->wavelength_list_gpu,
      spectrum_bands,
      use_filter_transmission);
  }

  cudaDeviceSynchronize(); 
  gpuErrchk( cudaPeekAtLastError() );
}



}


