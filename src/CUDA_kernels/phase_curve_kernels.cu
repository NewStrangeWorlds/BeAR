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
