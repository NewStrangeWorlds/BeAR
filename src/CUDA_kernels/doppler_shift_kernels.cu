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


#include "doppler_shift_kernels.h"
#include "error_check.h"

#include <cmath>


namespace bear {


// Doppler shift via sub-pixel interpolation on a log-lambda grid.
// Each thread handles one output pixel.
__global__
void dopplerShiftKernel(
  const float* __restrict__ spectrum_in,
  float*       __restrict__ spectrum_out,
  const int n_pixels,
  const float shift_pixels)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n_pixels) return;

  // Source position (fractional pixel in input)
  const float src = (float)i - shift_pixels;

  // Boundary: zero outside array
  if (src < 0.0f || src >= (float)(n_pixels - 1))
  {
    spectrum_out[i] = 0.0f;
    return;
  }

  const int idx = (int)src;
  const float frac = src - (float)idx;

  spectrum_out[i] = (1.0f - frac) * spectrum_in[idx]
                   + frac * spectrum_in[idx + 1];
}


__host__
void dopplerShiftGPU(
  const float* spectrum_in_dev,
  float* spectrum_out_dev,
  int n_pixels,
  float shift_pixels)
{
  const int threads = 256;
  const int blocks = (n_pixels + threads - 1) / threads;

  dopplerShiftKernel<<<blocks, threads>>>(
    spectrum_in_dev, spectrum_out_dev, n_pixels, shift_pixels);

  gpuErrchk(cudaDeviceSynchronize());
  gpuErrchk(cudaPeekAtLastError());
}


// Interpolate a model spectrum onto a target wavelength grid.
// Uses binary search + linear interpolation per target point.
__global__
void interpolateSpectrumKernel(
  const float* __restrict__ model_spectrum,
  const float* __restrict__ model_wavelengths,
  const int nb_model_points,
  const float* __restrict__ target_wavelengths,
  float* __restrict__ target_spectrum,
  const int nb_target_points)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nb_target_points) return;

  const float wl = target_wavelengths[i];

  // Binary search for the interval containing wl
  int lo = 0, hi = nb_model_points - 1;

  if (wl <= model_wavelengths[lo] || wl >= model_wavelengths[hi])
  {
    target_spectrum[i] = 0.0f;
    return;
  }

  while (hi - lo > 1)
  {
    int mid = (lo + hi) / 2;
    if (model_wavelengths[mid] <= wl)
      lo = mid;
    else
      hi = mid;
  }

  float t = (wl - model_wavelengths[lo])
           / (model_wavelengths[hi] - model_wavelengths[lo]);
  target_spectrum[i] = (1.0f - t) * model_spectrum[lo]
                      + t * model_spectrum[hi];
}


__host__
void interpolateSpectrumGPU(
  const float* model_spectrum_dev,
  const float* model_wavelengths_dev,
  int nb_model_points,
  const float* target_wavelengths_dev,
  float* target_spectrum_dev,
  int nb_target_points)
{
  const int threads = 256;
  const int blocks = (nb_target_points + threads - 1) / threads;

  interpolateSpectrumKernel<<<blocks, threads>>>(
    model_spectrum_dev, model_wavelengths_dev, nb_model_points,
    target_wavelengths_dev, target_spectrum_dev, nb_target_points);

  gpuErrchk(cudaDeviceSynchronize());
  gpuErrchk(cudaPeekAtLastError());
}


} // namespace bear
