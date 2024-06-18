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
#include <cmath>
#include <stdio.h>

#include "error_check.h"


namespace helios{


__global__ void addShiftToSpectrumDev(
  double* spectrum_bands,
  const unsigned int start_index,
  const double spectrum_shift,
  const int nb_points)
{

  //the thread index tid is the observation point
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_points; tid += blockDim.x * gridDim.x)
    spectrum_bands[start_index + tid] += spectrum_shift;

}



//converts layer optical depths to level-based extinction coefficients
__host__ void Observation::addShiftToSpectrumGPU(
  double* spectrum_bands,
  const unsigned int start_index,
  const double spectrum_shift)
{
  int threads = 256;
  int nb_points = nbPoints();

  int blocks = nb_points / threads;
  if (nb_points % threads) blocks++;

  addShiftToSpectrumDev<<<blocks,threads>>>(
    spectrum_bands,
    start_index,
    spectrum_shift,
    nb_points);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}


}
