// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

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
  float* stellar_spectrum,
  float* spectrum)
{
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_wavenumbers; tid += blockDim.x * gridDim.x)
  {
    spectrum[tid] = stellar_spectrum[tid];
  }
}



__host__ void StarSpectrumFile::calcFluxGPU(
  const std::vector<double>& parameter,
  float* spectrum_gpu)
{
  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  if (spectrum_dev == nullptr)
    moveToDevice(spectrum_dev, spectrum);

  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;
  
  const double effective_temperature = parameter[0];
  
  copyStellarSpectrum<<<blocks,threads>>>(
    nb_spectral_points,
    spectrum_dev,
    spectrum_gpu);

  CUDA_CHECK_AFTER_KERNEL();
}



}
 
