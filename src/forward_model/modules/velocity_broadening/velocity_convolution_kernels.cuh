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
 * velocity_convolution_kernels.cuh
 *
 * GPU velocity-space convolution for the high-resolution broadening modules
 * (velocity_broadening, phase_resolved_broadening).  Header-only, because the
 * templated kernel is instantiated in both modules' .cu files.
 *
 * The high-resolution grid is a subset of the opacity wavenumber points, so its
 * velocity step is NOT constant (at R = 2.5e5 on the 0.01 cm^-1 HELIOS-k grid it
 * varies by a factor of ~3 across 1.45-2.45 µm).  The convolution therefore works
 * with the actual step between neighbouring points, dv[i] = c |ln(λ_i/λ_{i+1})|:
 *
 *   out_i = Σ_j K(u_ij) w_j in_j / Σ_j K(u_ij) w_j,     |u_ij| <= half_width,
 *
 * with u_ij the velocity offset of point j from point i and w_j its trapezoidal
 * weight.  The kernel shape K only needs to be known up to a constant, since the
 * sum is normalised per output point.  This is the same scheme as the CPU code in
 * ../velocity_space_convolution.h.
 */


#ifndef _velocity_convolution_kernels_cuh
#define _velocity_convolution_kernels_cuh

#pragma once

#include <cmath>

#include "../../../CUDA_kernels/error_check.h"


namespace bear {


// Threads per block for the one-thread-per-output-point convolution kernels.
static constexpr int HIGHRES_BLOCK_SIZE = 128;


// Rotational broadening profile of Gray (2005, Eq. 18.14), unnormalised:
//   K(u) = 2(1-ε) sqrt(1-x²) + (π/2) ε (1-x²),  x = u / vsini,  |x| < 1
struct RotationalProfileGPU {
  float inv_vsini;
  float c1;
  float c2;

  __device__ __forceinline__ float operator()(const float u) const {
    const float x = u * inv_vsini;
    if (fabsf(x) >= 1.0f) return 0.0f;
    const float q = 1.0f - x * x;
    return c1 * sqrtf(q) + c2 * q;
  }
};


// Gaussian profile exp(-u² / 2σ²), unnormalised.
struct GaussianProfileGPU {
  float inv_2sig2;

  __device__ __forceinline__ float operator()(const float u) const {
    return expf(-u * u * inv_2sig2);
  }
};


// Profile tabulated at velocities v_min + k step (k = 0..size-1), linearly
// interpolated, zero outside the table.
struct TabulatedProfileGPU {
  const float* table;
  int size;
  float v_min;
  float inv_step;

  __device__ __forceinline__ float operator()(const float u) const {
    const float x = (u - v_min) * inv_step;
    if (x < 0.0f) return 0.0f;
    const int k = (int)x;
    if (k + 1 >= size) return 0.0f;
    const float t = x - (float)k;
    return (1.0f - t) * table[k] + t * table[k + 1];
  }
};


// One thread per output point.  The windows are a few tens of points wide, so a
// serial loop per thread is efficient, and adjacent threads read adjacent windows.
template <typename Profile>
__global__
void convolveVelocitySpaceKernel(
  const float* __restrict__ spectrum_in,
  float*       __restrict__ spectrum_out,
  const int                 n,
  const float* __restrict__ dv,
  const float               half_width,
  const Profile             profile)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= n) return;

  const float w_i = 0.5f * ((i > 0 ? dv[i - 1] : 0.0f) + (i < n - 1 ? dv[i] : 0.0f));

  float k = profile(0.0f) * w_i;
  float sum = k * spectrum_in[i];
  float norm = k;

  float u = 0.0f;

  for (int j = i + 1; j < n; ++j)
  {
    u += dv[j - 1];
    if (u > half_width) break;

    const float w_j = 0.5f * (dv[j - 1] + (j < n - 1 ? dv[j] : 0.0f));
    k = profile(u) * w_j;
    sum += k * spectrum_in[j];
    norm += k;
  }

  u = 0.0f;

  for (int j = i - 1; j >= 0; --j)
  {
    u += dv[j];
    if (u > half_width) break;

    const float w_j = 0.5f * ((j > 0 ? dv[j - 1] : 0.0f) + dv[j]);
    k = profile(-u) * w_j;
    sum += k * spectrum_in[j];
    norm += k;
  }

  spectrum_out[i] = norm > 0.0f ? sum / norm : spectrum_in[i];
}


template <typename Profile>
inline void launchVelocitySpaceConvolution(
  const float* spectrum_in_dev,
  float*       spectrum_out_dev,
  const int    n,
  const float* dv_dev,
  const float  half_width,
  const Profile& profile)
{
  const int blocks = (n + HIGHRES_BLOCK_SIZE - 1) / HIGHRES_BLOCK_SIZE;

  convolveVelocitySpaceKernel<<<blocks, HIGHRES_BLOCK_SIZE>>>(
    spectrum_in_dev, spectrum_out_dev, n, dv_dev, half_width, profile);

  CUDA_CHECK_AFTER_KERNEL();
}


} // namespace bear


#endif
