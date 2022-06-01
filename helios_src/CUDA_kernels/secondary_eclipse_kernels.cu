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


#include "../forward_model/secondary_eclipse/secondary_eclipse.h"

#include <iostream>
#include <cmath>
#include <stdio.h>

#include "error_check.h"
#include "reduce_kernels.h"
#include "../additional/physical_const.h"


namespace helios{


//computes the log likelihood
//each thread calculates one data point, the results are then summed by each block via a blockReduce method and finally all collected by thread 0 
__global__ void secondaryEclipseDevice(double* secondary_eclipse, double* planet_spectrum, double* stellar_spectrum,
                                       const int nb_points, const double radius_ratio, double* albedo_contribution)
{
  
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < nb_points; i += blockDim.x * gridDim.x)
  {

    secondary_eclipse[i] = planet_spectrum[i]/stellar_spectrum[i] * radius_ratio*radius_ratio * 1e6 + albedo_contribution[i]*1e6;

    //printf("%d  %1.5e  %1.5e  %1.5e  %f\n",  i, secondary_eclipse[i], planet_spectrum[i], stellar_spectrum[i], radius_ratio);
  
  }

}





__host__ void SecondaryEclipseModel::calcSecondaryEclipseGPU(double* secondary_eclipse, double* planet_spectrum_bands, double* stellar_spectrum_bands,
                                                             const int nb_points, const double radius_ratio, double* albedo_contribution)
{

  int threads = 256;

  int blocks = nb_points / threads;
  if (nb_points % threads) blocks++;


  secondaryEclipseDevice<<<blocks,threads>>>(secondary_eclipse, planet_spectrum_bands, stellar_spectrum_bands, nb_points,
                                             radius_ratio, albedo_contribution);

  
  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );

  //exit(0);
}


}
