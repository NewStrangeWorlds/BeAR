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


#include "../observations/observations.h"

#include <iostream>
#include <vector>
#include "math.h"
#include <stdio.h>
#include <new>

#include "../additional/physical_const.h"
#include "error_check.h"
#include "reduce_kernels.h"


namespace bear{


__global__ void applyFilterResponseDevice(
  const double* wavenumber, 
  double* spectrum, 
  double* filter_response_function, 
  double* filter_response_weight, 
  const double filter_normalisation,
  const int nb_points, 
  double* spectrum_filter)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < nb_points; i += blockDim.x * gridDim.x)
  {
    spectrum_filter[i] = spectrum[i] * filter_response_function[i] * filter_response_weight[i] / filter_normalisation;
  }

}



__host__ void Observation::applyFilterResponseGPU(double* spectrum)
{
  int nb_spectral_points = spectral_grid->nbSpectralPoints();

  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  applyFilterResponseDevice<<<blocks,threads>>>(
    spectral_grid->wavenumber_list_gpu,
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
