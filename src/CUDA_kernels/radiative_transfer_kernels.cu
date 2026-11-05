// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <iostream>
#include <cmath>
#include <stdio.h>

#include "../radiative_transfer/radiative_transfer.h"

#include "error_check.h"


namespace bear{


__global__ void changeSpectrumUnitsDevice(
  float* spectrum,
  double* wavelengths,
  int nb_points)
{
  
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < nb_points; i += blockDim.x * gridDim.x)
  {
    spectrum[i] = static_cast<float>(static_cast<double>(spectrum[i])/wavelengths[i]/wavelengths[i]*10000.0);
  }
}



__host__ void RadiativeTransfer::changeSpectrumUnitsGPU(
  float* spectrum)
{
  int threads = 256;
  
  int nb_points = spectral_grid->nbSpectralPoints();
  int blocks = nb_points / threads;
  if (nb_points % threads) blocks++;

  changeSpectrumUnitsDevice<<<blocks,threads>>>(
    spectrum,
    spectral_grid->wavelength_list_gpu,
    nb_points);

  CUDA_CHECK_AFTER_KERNEL();
}


}
