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
 * convolution_kernels_phase_resolved.cu
 *
 * GPU convolution with a pre-computed 1D broadening kernel for the
 * phase-resolved broadening module (Brogi et al. 2016).
 */


#include <cmath>
#include <cstdio>

#include "../../../CUDA_kernels/error_check.h"
#include "../../../CUDA_kernels/reduce_kernels.h"
#include "phase_resolved_convolution.h"


namespace bear {


static constexpr int PHASE_BLOCK_SIZE = 128;


// One block per output pixel.  Each thread handles a subset of the
// kernel window and the block reduces to the final sum.
__global__
void convolveWithPrecomputedKernel(
  const float* __restrict__ spectrum_in,
  float*       __restrict__ spectrum_out,
  const float* __restrict__ kernel_data,
  const int                 kernel_hw,
  const int                 n_pixels)
{
  const int i   = blockIdx.x;
  const int tid = threadIdx.x;

  if (i >= n_pixels) return;

  const int j_start = max(0,           i - kernel_hw);
  const int j_end   = min(n_pixels - 1, i + kernel_hw);

  float local_sum = 0.0f;

  for (int j = j_start + tid; j <= j_end; j += blockDim.x)
    local_sum += kernel_data[j - i + kernel_hw] * spectrum_in[j];

  local_sum = blockReduceSum(local_sum);

  if (tid == 0)
    spectrum_out[i] = local_sum;
}


__host__
void applyPrecomputedConvolutionGPU(
  const float* spectrum_in_dev,
  float*       spectrum_out_dev,
  const float* kernel_dev,
  int          kernel_hw,
  int          n_pixels)
{
  cudaGetLastError();

  convolveWithPrecomputedKernel<<<n_pixels, PHASE_BLOCK_SIZE>>>(
    spectrum_in_dev, spectrum_out_dev, kernel_dev, kernel_hw, n_pixels);

  CUDA_CHECK_AFTER_KERNEL();
}


}
