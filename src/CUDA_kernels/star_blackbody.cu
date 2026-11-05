// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <iostream>
#include <cmath>
#include <stdio.h>

#include "../forward_model/stellar_spectrum/star_blackbody.h"

#include "error_check.h"
#include "planck_function.h"
#include "../additional/physical_const.h"


namespace bear{


__global__ void starBlackBodyFlux(
  double* wavenumbers,
  int nb_wavenumbers,
  const double effective_temperature,
  float* flux)
{
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_wavenumbers; tid += blockDim.x * gridDim.x)
  {
    flux[tid] = static_cast<float>(planckFunction(
      effective_temperature,
      wavenumbers[tid]) * constants::pi * 1e-3);
  }
}



__host__ void StarBlackBody::calcFluxGPU(
  const std::vector<double>& parameter,
  float* spectrum_gpu)
{
  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();
  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  const double effective_temperature = parameter[0];


  starBlackBodyFlux<<<blocks,threads>>>(
    spectral_grid->wavenumber_list_gpu,
    nb_spectral_points,
    effective_temperature,
    spectrum_gpu);

  CUDA_CHECK_AFTER_KERNEL();
}


}
 
