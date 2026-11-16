// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <iostream>
#include <vector>
#include "math.h"
#include <stdio.h>
#include <new>

#include "../spectral_grid/spectral_band.h"
#include "../additional/physical_const.h"


#include "error_check.h"
#include <algorithm>
#include "reduce_kernels.h"
#include "../spectral_grid/spectral_grid.h"


namespace bear{


namespace {

//block size of the convolution kernels; the tile sizes computed on the host
//(SpectralBands::initDeviceMemory) refer to this value
constexpr int convolution_block_size = 256;

//largest dynamic shared memory the tiled kernel asks for (bytes)
constexpr int convolution_max_shared_bytes = 96 * 1024;


__device__ __forceinline__
float normalFactorFl(float sigma)
{
  return rsqrtf(2.0f * constants::pi) / sigma;  // 1/(σ√2π)
}


//convolution of the spectrum at the high-res point i with its Gaussian instrument
//profile over the grid points [start, end] (composite trapezoidal rule: each grid
//point is evaluated exactly once and weighted by the width of its adjacent intervals,
//the edge points only receive the width of their single interior interval)
//the wavelengths and the spectrum are read through accessors so that the same
//arithmetic serves the tiled (shared memory) and the direct (global memory) kernel
template <class Wavelength, class Spectrum>
__device__ __forceinline__
float convolvePoint(
  const int start,
  const int end,
  const float mu,
  const float sigma,
  const Wavelength& wavelength,
  const Spectrum& spectrum_value)
{
  const float norm = normalFactorFl(sigma);
  const float inv2sig2 = 1.0f / (2.0f * sigma * sigma);

  float sum = 0.0f;

  for (int j = start; j <= end; ++j)
  {
    const float wl = wavelength(j);
    const float d = wl - mu;
    const float g = __expf(-d * d * inv2sig2);

    float weight = 0.0f;
    if (j < end)   weight += wavelength(j + 1) - wl;
    if (j > start) weight += wl - wavelength(j - 1);

    sum += spectrum_value(j) * g * weight;
  }

  return fabsf(0.5f * norm * sum);
}

}



//one thread per high-res point of the observational range; the union of the
//convolution windows of the block is staged in shared memory (wavelengths as float
//and spectrum), so that the inner loops run from on-chip memory
__global__ 
void convolveSpectrumTiledDeviceFl(
  const float* __restrict__ spectrum,
  const int index_start,
  const int nb_points,
  const double* __restrict__ wavelengths,
  const double* __restrict__ band_sigma,
  const int* __restrict__ start_index,
  const int* __restrict__ end_index,
  const int tile_capacity,
  float* __restrict__ convolved_spectrum)
{
  extern __shared__ float shared_data[];
  float* s_wavelength = shared_data;
  float* s_spectrum = shared_data + tile_capacity;

  __shared__ int s_tile_start;
  __shared__ int s_tile_end;

  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const bool active = i < nb_points;

  const int start = active ? start_index[i] : 0x7fffffff;
  const int end = active ? end_index[i] : -1;

  if (threadIdx.x == 0)
  {
    s_tile_start = 0x7fffffff;
    s_tile_end = -1;
  }

  __syncthreads();

  atomicMin(&s_tile_start, start);
  atomicMax(&s_tile_end, end);

  __syncthreads();

  const int tile_start = s_tile_start;
  const int tile_size = s_tile_end - tile_start + 1;

  for (int k = threadIdx.x; k < tile_size; k += blockDim.x)
  {
    s_wavelength[k] = static_cast<float>(wavelengths[tile_start + k]);
    s_spectrum[k] = spectrum[tile_start + k];
  }

  __syncthreads();

  if (!active) return;

  const float sigma = static_cast<float>(band_sigma[i + index_start]);

  if (start == end || sigma == 0.0f)
  {
    convolved_spectrum[i + index_start] = spectrum[i + index_start];
    return;
  }

  const float mu = static_cast<float>(wavelengths[i + index_start]);

  convolved_spectrum[i + index_start] = convolvePoint(
    start, end, mu, sigma,
    [&](const int j) { return s_wavelength[j - tile_start]; },
    [&](const int j) { return s_spectrum[j - tile_start]; });
}



//one thread per high-res point, windows read from global memory; used when the
//windows are too wide for the shared-memory tile
__global__ 
void convolveSpectrumDirectDeviceFl(
  const float* __restrict__ spectrum,
  const int index_start,
  const int nb_points,
  const double* __restrict__ wavelengths,
  const double* __restrict__ band_sigma,
  const int* __restrict__ start_index,
  const int* __restrict__ end_index,
  float* __restrict__ convolved_spectrum)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= nb_points) return;

  const int start = start_index[i];
  const int end = end_index[i];
  const float sigma = static_cast<float>(band_sigma[i + index_start]);

  if (start == end || sigma == 0.0f)
  {
    convolved_spectrum[i + index_start] = spectrum[i + index_start];
    return;
  }

  const float mu = static_cast<float>(wavelengths[i + index_start]);

  convolved_spectrum[i + index_start] = convolvePoint(
    start, end, mu, sigma,
    [&](const int j) { return static_cast<float>(wavelengths[j]); },
    [&](const int j) { return spectrum[j]; });
}



//largest number of grid points any block of the tiled kernel has to stage,
//evaluated on the host from the convolution windows
int SpectralBands::convolutionTileSize(
  const std::vector<int>& convolution_start,
  const std::vector<int>& convolution_end)
{
  const int nb_points = convolution_start.size();
  int max_tile = 0;

  for (int block_start = 0; block_start < nb_points; block_start += convolution_block_size)
  {
    const int block_end = std::min(block_start + convolution_block_size, nb_points);

    int tile_start = convolution_start[block_start];
    int tile_end = convolution_end[block_start];

    for (int i = block_start; i < block_end; ++i)
    {
      tile_start = std::min(tile_start, convolution_start[i]);
      tile_end = std::max(tile_end, convolution_end[i]);
    }

    max_tile = std::max(max_tile, tile_end - tile_start + 1);
  }

  return max_tile;
}



__host__ 
void SpectralBands::convolveSpectrumGPU(
  float* spectrum,
  float* spectrum_processed_dev)
{
  const int nb_high_res_points = 
    obs_index_range.second - obs_index_range.first + 1;

  const int threads = convolution_block_size;
  const int blocks = (nb_high_res_points + threads - 1) / threads;

  const size_t shared_bytes = 2 * static_cast<size_t>(convolution_max_tile) * sizeof(float);

  if (shared_bytes <= convolution_max_shared_bytes)
  {
    static size_t attribute_bytes = 0;

    if (shared_bytes > 48 * 1024 && shared_bytes > attribute_bytes)
    {
      gpuErrchk(cudaFuncSetAttribute(
        convolveSpectrumTiledDeviceFl,
        cudaFuncAttributeMaxDynamicSharedMemorySize,
        convolution_max_shared_bytes));
      attribute_bytes = convolution_max_shared_bytes;
    }

    convolveSpectrumTiledDeviceFl<<<blocks, threads, shared_bytes>>>(
      spectrum, 
      obs_index_range.first,
      nb_high_res_points,
      spectral_grid->wavelength_list_gpu, 
      instrument_profile_sigma_dev,
      convolution_start_dev, 
      convolution_end_dev,
      convolution_max_tile,
      spectrum_processed_dev);
  }
  else
  {
    convolveSpectrumDirectDeviceFl<<<blocks, threads>>>(
      spectrum, 
      obs_index_range.first,
      nb_high_res_points,
      spectral_grid->wavelength_list_gpu, 
      instrument_profile_sigma_dev,
      convolution_start_dev, 
      convolution_end_dev,
      spectrum_processed_dev);
  }

  CUDA_CHECK_AFTER_KERNEL();
}



}
