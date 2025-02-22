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
#include <cmath>
#include <stdio.h>

#include "../forward_model/secondary_eclipse/secondary_eclipse.h"

#include "error_check.h"
#include "reduce_kernels.h"
#include "../additional/physical_const.h"


namespace bear{


__global__ void OccultationDevice(
  double* secondary_eclipse,
  double* planet_spectrum,
  const double* stellar_spectrum,
  const int nb_points,
  const double radius_ratio, 
  const double* albedo_contribution)
{
  
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < nb_points; i += blockDim.x * gridDim.x)
  {
    secondary_eclipse[i] = planet_spectrum[i]/stellar_spectrum[i] * radius_ratio*radius_ratio * 1e6; 

    if (albedo_contribution != nullptr)
      secondary_eclipse[i] += albedo_contribution[i]*1e6;
  }
}



__host__ void OccultationModel::calcOccultationGPU(
  double* secondary_eclipse,
  double* planet_spectrum,
  const double* stellar_spectrum,
  const int nb_points,
  const double radius_ratio,
  const double* albedo_contribution)
{
  int threads = 256;

  int blocks = nb_points / threads;
  if (nb_points % threads) blocks++;


  OccultationDevice<<<blocks,threads>>>(
    secondary_eclipse,
    planet_spectrum,
    stellar_spectrum,
    nb_points,
    radius_ratio,
    albedo_contribution);


  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}


}
