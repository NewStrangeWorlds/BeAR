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


/*
 * convolution_kernels_highres.cu
 *
 * GPU spectral broadening for high-resolution spectroscopy.
 * Two-pass pipeline: rotational broadening (Gray 2005) then
 * instrumental Gaussian broadening, both in velocity space
 * on a constant-resolution (log-lambda) grid.
 */


#include <cmath>
#include <stdio.h>

#include "../../../CUDA_kernels/error_check.h"
#include "../../../CUDA_kernels/reduce_kernels.h"
#include "../../../CUDA_kernels/data_management_kernels.h"
#include "velocity_convolution_kernels.cuh"
#include "../../../CUDA_kernels/highres_convolution.h"


namespace bear {


// One thread per output pixel.  The convolution windows are only a few to a
// few tens of pixels wide (5 sigma of the instrumental profile, or vsini, in
// pixel units), so a serial loop per thread is far more efficient than a
// block-wide reduction: no idle threads, no reduction overhead, and adjacent
// threads read adjacent windows so the loads stay coalesced.
__global__
void convolveGaussianHRKernel(
  const float* __restrict__ spectrum_in,
  float*       __restrict__ spectrum_out,
  const int                 n_pixels,
  const float               sigma_pixels,
  const int                 half_width)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= n_pixels) return;

  const float inv_2sig2 = 1.0f / (2.0f * sigma_pixels * sigma_pixels);
  const float norm      = rsqrtf(2.0f * (float)M_PI) / sigma_pixels;

  const int j_start = max(0,           i - half_width);
  const int j_end   = min(n_pixels - 1, i + half_width);

  float local_sum = 0.0f;

  for (int j = j_start; j <= j_end; ++j)
  {
    float dx = (float)(j - i);
    local_sum += gaussianKernelHR(dx, inv_2sig2, norm) * spectrum_in[j];
  }

  spectrum_out[i] = local_sum;
}


__global__
void convolveRotationalHRKernel(
  const float* __restrict__ spectrum_in,
  float*       __restrict__ spectrum_out,
  const int                 n_pixels,
  const float               vsini_pixels,
  const float               epsilon)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= n_pixels) return;

  const int half_width = (int)ceilf(vsini_pixels);
  const int j_start    = max(0,           i - half_width);
  const int j_end      = min(n_pixels - 1, i + half_width);

  float local_sum = 0.0f;

  for (int j = j_start; j <= j_end; ++j)
  {
    float dx = (float)(j - i);
    local_sum += rotationalKernelHR(dx, vsini_pixels, epsilon) * spectrum_in[j];
  }

  spectrum_out[i] = local_sum;
}


__host__
void applyHighResConvolutionGPU(
  float*  spectrum_in_dev,
  float*  spectrum_out_dev,
  int     n_pixels,
  double  sigma_kms,
  double  vsini_kms,
  double  delta_v_kms,
  double  epsilon,
  float*  temp_dev)
{
  cudaGetLastError();

  const int threads = HIGHRES_BLOCK_SIZE;
  const int blocks  = (n_pixels + threads - 1) / threads;

  const float sigma_pixels = (float)(sigma_kms / delta_v_kms);
  const int   half_width   = (int)ceil(5.0 * sigma_kms / delta_v_kms);

  const bool do_rotation = vsini_kms > 0.5 * delta_v_kms;
  const bool do_gaussian = sigma_pixels > 0.01f;

  bool allocated_temp = false;

  if (do_rotation && do_gaussian)
  {
    if (temp_dev == nullptr)
    {
      allocateOnDevice(temp_dev, (size_t)n_pixels);
      allocated_temp = true;
    }

    const float vsini_pixels = (float)(vsini_kms / delta_v_kms);

    convolveRotationalHRKernel<<<blocks, threads>>>(
      spectrum_in_dev, temp_dev, n_pixels, vsini_pixels, (float)epsilon);

    CUDA_CHECK_AFTER_KERNEL();

    convolveGaussianHRKernel<<<blocks, threads>>>(
      temp_dev, spectrum_out_dev, n_pixels, sigma_pixels, half_width);
  }
  else if (do_rotation)
  {
    const float vsini_pixels = (float)(vsini_kms / delta_v_kms);

    convolveRotationalHRKernel<<<blocks, threads>>>(
      spectrum_in_dev, spectrum_out_dev, n_pixels, vsini_pixels, (float)epsilon);
  }
  else if (do_gaussian)
  {
    convolveGaussianHRKernel<<<blocks, threads>>>(
      spectrum_in_dev, spectrum_out_dev, n_pixels, sigma_pixels, half_width);
  }
  else
  {
    // No broadening — copy input to output
    copyOnDevice(spectrum_out_dev, spectrum_in_dev, (size_t)n_pixels);
  }

  CUDA_CHECK_AFTER_KERNEL();

  if (allocated_temp)
    deleteFromDevice(temp_dev);
}


} // namespace bear
