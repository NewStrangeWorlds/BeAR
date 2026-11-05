// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include "../observations/observations.h"


#include <iostream>
#include <cmath>
#include <stdio.h>

#include "error_check.h"


namespace bear{


__global__ 
void addShiftToSpectrumDev(
  float* __restrict__ spectrum_bands,
  const double spectrum_shift,
  const int nb_points)
{

  //the thread index tid is the observation point
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_points; tid += blockDim.x * gridDim.x)
    spectrum_bands[tid] += static_cast<float>(spectrum_shift);

}


__host__ 
void Observation::addShiftToSpectrumGPU(
  float* spectrum_obs,
  const double spectrum_shift)
{
  int threads = 256;
  int nb_points = nbPoints();

  int blocks = nb_points / threads;
  if (nb_points % threads) blocks++;

  addShiftToSpectrumDev<<<blocks,threads>>>(
    spectrum_obs,
    spectrum_shift,
    nb_points);

  CUDA_CHECK_AFTER_KERNEL();
}


}
