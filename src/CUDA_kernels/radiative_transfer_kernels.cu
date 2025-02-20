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

#include "../radiative_transfer/radiative_transfer.h"

#include "error_check.h"


namespace bear{


__global__ void changeSpectrumUnitsDevice(
  double* spectrum,
  double* wavelengths,
  int nb_points)
{
  
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < nb_points; i += blockDim.x * gridDim.x)
  {
    spectrum[i] = spectrum[i]/wavelengths[i]/wavelengths[i]*10000.0;
  }
}



__host__ void RadiativeTransfer::changeSpectrumUnitsGPU(
  double* spectrum)
{
  int threads = 256;
  
  int nb_points = spectral_grid->nbSpectralPoints();
  int blocks = nb_points / threads;
  if (nb_points % threads) blocks++;

  changeSpectrumUnitsDevice<<<blocks,threads>>>(
    spectrum,
    spectral_grid->wavelength_list_gpu,
    nb_points);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}


}
