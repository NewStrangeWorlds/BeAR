// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <iostream>
#include <cmath>
#include <stdio.h>

#include "../forward_model/secondary_eclipse/secondary_eclipse.h"

#include "error_check.h"
#include "reduce_kernels.h"
#include "../additional/physical_const.h"


namespace bear{


__global__ void OccultationDevice(
  float* secondary_eclipse,
  float* planet_spectrum,
  const float* stellar_spectrum,
  const int nb_points,
  const double radius_ratio,
  const float* albedo_contribution)
{

  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < nb_points; i += blockDim.x * gridDim.x)
  {
    secondary_eclipse[i] = static_cast<float>(static_cast<double>(planet_spectrum[i])/static_cast<double>(stellar_spectrum[i]) * radius_ratio*radius_ratio * 1e6);

    if (albedo_contribution != nullptr)
      secondary_eclipse[i] += static_cast<float>(static_cast<double>(albedo_contribution[i])*1e6);
  }
}



__host__ void OccultationModel::calcOccultationGPU(
  float* secondary_eclipse,
  float* planet_spectrum,
  const float* stellar_spectrum,
  const int nb_points,
  const double radius_ratio,
  const float* albedo_contribution)
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


  CUDA_CHECK_AFTER_KERNEL();
}


}
