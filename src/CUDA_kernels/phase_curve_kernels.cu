// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include "../forward_model/phase_curve/phase_curve.h"
#include "error_check.h"
#include "data_management_kernels.h"


namespace bear{


// Divide planet spectrum by stellar spectrum and scale by (Rp/Rs)^2
// to produce the dimensionless Fp/Fs ratio in-place.
__global__ void normaliseFpFsDevice(
  float*       planet_spectrum,
  const float* stellar_spectrum,
  const int    nb_points,
  const float  radius_ratio_squared)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < nb_points;
       i += blockDim.x * gridDim.x)
  {
    planet_spectrum[i] = static_cast<float>(
      static_cast<double>(planet_spectrum[i])
      / static_cast<double>(stellar_spectrum[i])
      * static_cast<double>(radius_ratio_squared));
  }
}


__global__ void calcOccultationLowResDevice(
  float*       output,
  const float* planet_spectrum,
  const float* stellar_spectrum,
  const int    nb_points,
  const float  radius_ratio_squared)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < nb_points;
       i += blockDim.x * gridDim.x)
  {
    output[i] = static_cast<float>(
      static_cast<double>(planet_spectrum[i])
      / static_cast<double>(stellar_spectrum[i])
      * static_cast<double>(radius_ratio_squared)
      * 1.0e6);
  }
}


__host__ void PhaseCurveModel::calcOccultationLowResGPU(
  float*       output,
  const float* planet_spectrum,
  const float* stellar_spectrum,
  const int    nb_points,
  const float  radius_ratio_squared)
{
  const int threads = 256;
  const int blocks  = (nb_points + threads - 1) / threads;

  calcOccultationLowResDevice<<<blocks, threads>>>(
    output,
    planet_spectrum,
    stellar_spectrum,
    nb_points,
    radius_ratio_squared);

  CUDA_CHECK_AFTER_KERNEL();
}


__host__ void PhaseCurveModel::normaliseFpFsGPU(
  float*       planet_spectrum,
  const float* stellar_spectrum,
  const int    nb_points,
  const float  radius_ratio_squared)
{
  const int threads = 256;
  const int blocks  = (nb_points + threads - 1) / threads;

  normaliseFpFsDevice<<<blocks, threads>>>(
    planet_spectrum,
    stellar_spectrum,
    nb_points,
    radius_ratio_squared);

  CUDA_CHECK_AFTER_KERNEL();
}


}
