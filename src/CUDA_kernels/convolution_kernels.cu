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
#include <vector>
#include "math.h"
#include <stdio.h>
#include <new>

#include "../spectral_grid/spectral_band.h"
#include "../additional/physical_const.h"


#include "error_check.h"
#include "reduce_kernels.h"
#include "../spectral_grid/spectral_grid.h"


namespace bear{


__device__ __forceinline__
float normalFactorFl(float sigma)
{
  return rsqrtf(2.0f * constants::pi) / sigma;  // 1/(σ√2π)
}


__global__ 
void convolveSpectrumDeviceFl(
  const float* __restrict__ spectrum,
  const int index_start,
  const double* __restrict__ wavelengths,
  const double* __restrict__ band_sigma,
  const int* __restrict__ start_index,
  const int* __restrict__ end_index,
  float* __restrict__ convolved_spectrum)
{
  const int i = blockIdx.x;
  const int tid = threadIdx.x;

  const float mu = wavelengths[i + index_start];
  const float sigma = band_sigma[i + index_start];

  const int start = start_index[i];
  const int end   = end_index[i];

  if (start == end || sigma == 0.0) 
  {
    if (tid == 0)
      convolved_spectrum[i + index_start] = spectrum[i + index_start];
    return;
  }

  const float norm = normalFactorFl(sigma);
  const float inv2sig2 = 1.0f / (2.0f * sigma * sigma);

  float local_sum = 0.0f;

  //composite trapezoidal rule: each grid point is evaluated exactly once and
  //weighted by the width of its adjacent intervals
  //(the edge points only receive the width of their single interior interval)
  for (int j = start + tid; j <= end; j += blockDim.x)
  {
    const float wl = wavelengths[j];
    const float d = wl - mu;
    const float g = __expf(-d * d * inv2sig2);

    float weight = 0.0f;
    if (j < end)   weight += (float)wavelengths[j + 1] - wl;
    if (j > start) weight += wl - (float)wavelengths[j - 1];

    local_sum += spectrum[j] * g * weight;
  }

  local_sum = blockReduceSum(local_sum);

  if (tid == 0)
    convolved_spectrum[i + index_start] = fabs(0.5f * norm * local_sum);
}



__host__ 
void SpectralBands::convolveSpectrumGPU(
  float* spectrum,
  float* spectrum_processed_dev)
{
  const size_t nb_high_res_points = 
    obs_index_range.second - obs_index_range.first + 1;

  int threads = 128;
  int blocks = nb_high_res_points;

  convolveSpectrumDeviceFl<<<blocks,threads>>>(
    spectrum, 
    obs_index_range.first,
    spectral_grid->wavelength_list_gpu, 
    instrument_profile_sigma_dev,
    convolution_start_dev, 
    convolution_end_dev,
    spectrum_processed_dev);

  CUDA_CHECK_AFTER_KERNEL();
}



}
