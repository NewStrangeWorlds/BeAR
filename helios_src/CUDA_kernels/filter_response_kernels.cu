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



#include <iostream>
#include <vector>
#include "math.h"
#include <stdio.h>
#include <new>

#include "../additional/physical_const.h"
#include "../observations/observations.h"
#include "../retrieval/retrieval.h"

#include "error_check.h"
#include "reduce_kernels.h"


namespace helios{



__global__ void applyFilterResponseDevice(
  const double* wavenumber, 
  double* spectrum, 
  double* filter_response_function, 
  double* filter_response_weight, 
  const double filter_normalisation,
  const int nb_points, 
  double* convolved_spectrum)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < nb_points; i += blockDim.x * gridDim.x)
  {
    convolved_spectrum[i] = spectrum[i] * filter_response_function[i] * filter_response_weight[i] / filter_normalisation;
  }

}



__host__ void Observation::applyFilterResponseGPU(double* spectrum)
{
  int nb_spectral_points = retrieval->spectral_grid.nbSpectralPoints();

  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  applyFilterResponseDevice<<<blocks,threads>>>(
    retrieval->spectral_grid.wavenumber_list_gpu,
    spectrum,
    filter_response_gpu,
    filter_response_weight_gpu,
    filter_response_normalisation,
    nb_spectral_points, 
    spectrum_filter_dev);


  cudaDeviceSynchronize(); 
  gpuErrchk( cudaPeekAtLastError() ); 
}



}
