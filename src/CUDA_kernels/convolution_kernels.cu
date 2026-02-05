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
  const double* __restrict__ spectrum,
  const int index_start,
  const double* __restrict__ wavelengths,
  const double* __restrict__ band_sigma,
  const int* __restrict__ start_index,
  const int* __restrict__ end_index,
  double* __restrict__ convolved_spectrum)
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
  const float inv2sig2 = 1.0 / (2.0 * sigma * sigma);

  float local_sum = 0.0;

  //each thread integrates part of the trapezoids directly
  for (int j = start + tid; j < end; j += blockDim.x)
  {
    float wl0 = wavelengths[j];
    float wl1 = wavelengths[j + 1];

    float d0 = wl0 - mu;
    float d1 = wl1 - mu;

    float g0 = norm * exp(-d0 * d0 * inv2sig2);
    float g1 = norm * exp(-d1 * d1 * inv2sig2);

    float s0 = spectrum[j] * g0;
    float s1 = spectrum[j + 1] * g1;

    local_sum += (s0 + s1) * (wl1 - wl0);
  }

  local_sum = blockReduceSum(local_sum);

  if (tid == 0)
    convolved_spectrum[i + index_start] = fabs(0.5 * local_sum);
}




__device__ __forceinline__
double normalFactor(double sigma)
{
  return rsqrt(2.0 * constants::pi) / sigma;  // 1/(σ√2π)
}


__global__ 
void convolveSpectrumDevice(
  const double* __restrict__ spectrum,
  const int index_start,
  const double* __restrict__ wavelengths,
  const double* __restrict__ band_sigma,
  const int* __restrict__ start_index,
  const int* __restrict__ end_index,
  double* __restrict__ convolved_spectrum)
{
  const int i = blockIdx.x;
  const int tid = threadIdx.x;

  const double mu = wavelengths[i + index_start];
  const double sigma = band_sigma[i + index_start];

  const int start = start_index[i];
  const int end   = end_index[i];

  if (start == end || sigma == 0.0) 
  {
    if (tid == 0)
      convolved_spectrum[i + index_start] = spectrum[i + index_start];
    return;
  }

  const float norm = normalFactor(sigma);
  const float inv2sig2 = 1.0 / (2.0 * sigma * sigma);

  double local_sum = 0.0;

  //each thread integrates part of the trapezoids directly
  for (int j = start + tid; j < end; j += blockDim.x)
  {
    double wl0 = wavelengths[j];
    double wl1 = wavelengths[j + 1];

    double d0 = wl0 - mu;
    double d1 = wl1 - mu;

    double g0 = norm * exp(-d0 * d0 * inv2sig2);
    double g1 = norm * exp(-d1 * d1 * inv2sig2);

    double s0 = spectrum[j]     * g0;
    double s1 = spectrum[j + 1] * g1;

    local_sum += (s0 + s1) * (wl1 - wl0);
  }

  local_sum = blockReduceSum(local_sum);

  if (tid == 0)
    convolved_spectrum[i + index_start] = fabs(0.5 * local_sum);
}



__host__ 
void SpectralBands::convolveSpectrumGPU(
  double* spectrum, 
  double* spectrum_processed_dev)
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

  cudaDeviceSynchronize(); 
  gpuErrchk( cudaPeekAtLastError() ); 
}



}
