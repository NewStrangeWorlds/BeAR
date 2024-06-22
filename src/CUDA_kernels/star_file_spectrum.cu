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

#include "../forward_model/stellar_spectrum/star_file_spectrum.h"

#include "error_check.h"
#include "planck_function.h"
#include "../additional/physical_const.h"


namespace bear{


__global__ void copyStellarSpectrum(
  int nb_wavenumbers,
  double* stellar_spectrum,
  double* spectrum)
{
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_wavenumbers; tid += blockDim.x * gridDim.x)
  {
    spectrum[tid] = stellar_spectrum[tid];
  }
}



__host__ void StarSpectrumFile::calcFluxGPU(
  const std::vector<double>& parameter,
  double* spectrum_gpu)
{
  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  if (spectrum_dev == nullptr)
    moveToDevice(spectrum_dev, spectrum, true);

  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;
  
  const double effective_temperature = parameter[0];
  
  copyStellarSpectrum<<<blocks,threads>>>(
    nb_spectral_points,
    spectrum_dev,
    spectrum_gpu);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}



}
 
