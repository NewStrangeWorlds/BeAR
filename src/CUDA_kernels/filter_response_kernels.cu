// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

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
  float* spectrum,
  double* filter_response_function,
  double* filter_response_weight,
  const double filter_normalisation,
  const int nb_points,
  float* spectrum_filter)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < nb_points; i += blockDim.x * gridDim.x)
  {
    spectrum_filter[i] = static_cast<float>(static_cast<double>(spectrum[i]) * filter_response_function[i] * filter_response_weight[i] / filter_normalisation);
  }

}



__host__ void Observation::applyFilterResponseGPU(float* spectrum)
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


  CUDA_CHECK_AFTER_KERNEL();
}


}
