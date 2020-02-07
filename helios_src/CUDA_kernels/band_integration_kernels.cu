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


#include "band_integration_kernels.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <new>

#include "data_management_kernels.h"
#include "error_check.h"
#include "reduce_kernels.h"



namespace helios{


//every block reduces one band
//spectrum_high_res is the pointer to the array on the GPU
//the integration is done by a simple piece-wise trapezoidal rule
//note that the units of the high-res spectrum are in W m-2 cm, while the mean band values are in W m-2 mu-1
__global__ void bandIntegrationDevice(double* spectrum_high_res, int* band_indices, int* band_indices_size, int* band_start_index,
                                      int nb_bands,
                                      double* wavenumbers, double* wavelengths,
                                      double* spectrum_bands)
{
  double band_sum = 0;

  //the current band index
  const int i = blockIdx.x;

  //indices to navigate through the high-res spectrum
  const int start_index = band_start_index[i];
  const int end_index = band_start_index[i] + band_indices_size[i] - 1;
  const int band_size = band_indices_size[i];


  for (int j = threadIdx.x; j < band_size-1; j += blockDim.x)
  {
    const int index1 = band_indices[j + start_index + 1];
    const int index2 = band_indices[j + start_index];
      
    const double sum = (spectrum_high_res[band_indices[j + start_index + 1] ] + spectrum_high_res[band_indices[j + start_index]])
                       * (wavenumbers[index1] - wavenumbers[index2]);

    band_sum += sum;
  }
  

  __syncthreads();


  band_sum = blockReduceSum(band_sum);


  if (threadIdx.x == 0)
    spectrum_bands[blockIdx.x] = band_sum * 0.5 / (wavelengths[band_indices[start_index]] - wavelengths[band_indices[end_index]]);
}



//every block reduces one band
//spectra_high_res is the pointer to the pointer to the data array on the GPU 
//(there may be multiple convolved high-res spectra and we need to pick the correct one)
//the integration is done by a simple piece-wise trapezoidal rule
//note that the units of the high-res spectrum are in W m-2 cm, while the mean band values are in W m-2 mu-1
__global__ void bandIntegrationDevice(double** spectra_high_res, int* band_indices, int* band_indices_size, int* band_start_index,
                                      int nb_bands,
                                      double* wavenumbers, double* wavelengths,
                                      double* spectrum_bands)
{
  double band_sum = 0;

  //the current band index
  const int i = blockIdx.x;


  //pick the correct spectrum to integrate
  double* spectrum_high_res = spectra_high_res[i];


  //indices to navigate through the high-res spectrum
  const int start_index = band_start_index[i];
  const int end_index = band_start_index[i] + band_indices_size[i] - 1;
  const int band_size = band_indices_size[i];


  for (int j = threadIdx.x; j < band_size-1; j += blockDim.x)
  {
    const int index1 = band_indices[j + start_index + 1];
    const int index2 = band_indices[j + start_index];

    const double sum = (spectrum_high_res[band_indices[j + start_index + 1] ] + spectrum_high_res[band_indices[j + start_index]])
                     * (wavenumbers[index1] - wavenumbers[index2]);

    band_sum += sum;
  }


  __syncthreads();


  band_sum = blockReduceSum(band_sum);


  if (threadIdx.x == 0)
    spectrum_bands[blockIdx.x] = band_sum * 0.5 / (wavelengths[band_indices[start_index]] - wavelengths[band_indices[end_index]]);
}



//host function for the band integration of the high-res spectrum
//spectra_high_res contains the GPU pointer to the high-res spectrum
__host__ void bandIntegrationGPU(double* spectrum_high_res, int* band_indices, int* band_indices_size, int* band_start_index,
                                 int nb_bands,
                                 double* wavenumbers, double* wavelengths,
                                 double* spectrum_bands)
{
  int threads = 128;
  int blocks = nb_bands;


  bandIntegrationDevice<<<blocks,threads>>>(spectrum_high_res, band_indices, band_indices_size, band_start_index,
                                            nb_bands,
                                            wavenumbers, wavelengths,
                                            spectrum_bands);

 
  cudaDeviceSynchronize(); 
  gpuErrchk( cudaPeekAtLastError() ); 
}



//host function for the band integration of the high-res spectrum
//spectra_high_res contains a vector of GPU pointers, each with a differently convolved spectrum 
__host__ void bandIntegrationGPU(std::vector<double*> spectra_high_res, int* band_indices, int* band_indices_size, int* band_start_index,
                                 int nb_bands,
                                 double* wavenumbers, double* wavelengths,
                                 double* spectrum_bands)
{
  int threads = 128;
  int blocks = nb_bands;


  double** spectra_high_res_dev;

  moveToDevice(spectra_high_res_dev, spectra_high_res);

 

  bandIntegrationDevice<<<blocks,threads>>>(spectra_high_res_dev, band_indices, band_indices_size, band_start_index,
                                            nb_bands,
                                            wavenumbers, wavelengths,
                                            spectrum_bands);

  cudaDeviceSynchronize(); 
  gpuErrchk( cudaPeekAtLastError() ); 
  
  
  deleteFromDevice(spectra_high_res_dev);
}




}


