// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include "../forward_model/transmission/transmission.h"
#include "../forward_model/atmosphere/atmosphere.h"

#include <iostream>
#include <cmath>
#include <stdio.h>

#include "error_check.h"
#include "reduce_kernels.h"
#include "data_management_kernels.h"
#include "../additional/physical_const.h"


namespace bear{


namespace {

//threads of a block of the transmission kernel: 32 wavelengths times
//transmission_block_heights tangent heights that are processed in turn
constexpr int transmission_block_wavelengths = 32;
constexpr int transmission_block_heights = 16;

}


//geometric path lengths of the limb rays: path_lengths[t * nb_grid_points + i] is the
//length of the ray with tangent height t inside the layer between levels i and i+1
//(one hemisphere), for i >= t; evaluated once per model evaluation
__global__
void transmissionPathLengthKernel(
  const int nb_grid_points,
  const float bottom_radius,
  const float* __restrict__ altitude_ns,
  float* __restrict__ path_lengths)
{
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= nb_grid_points * nb_grid_points) return;

  const int t = idx / nb_grid_points;
  const int i = idx - t * nb_grid_points;

  if (i < t || i >= nb_grid_points - 1) return;

  const float b = bottom_radius + altitude_ns[t];
  const float a1 = bottom_radius + altitude_ns[i];
  const float a2 = bottom_radius + altitude_ns[i+1];

  float path = sqrtf(a2*a2 - b*b);

  if (i != t)
    path -= sqrtf(a1*a1 - b*b);

  path_lengths[idx] = path;
}



//transit radius per wavelength: a block handles 32 consecutive wavelengths
//the block first collects the total extinction coefficients (gas absorption, 
//scattering, and cloud extinction) of its wavelengths in shared memory and adds 
//those of adjacent levels, so that the tangent optical depths are evaluated without 
//any further access to the large coefficient arrays
//the cloud extinction coefficient of a level is the optical depth of the layer below
//divided by the layer's thickness
//scattering_coeff and cloud_optical_depth can be a nullptr; the former if the 
//absorption coefficients already contain the gas scattering
//its thread rows then take turns over the tangent heights and evaluate the tangent 
//optical depth and transmission of their rays, the first row finally integrates 
//the effective tangent height from the top down
//shared memory: (2*nb_grid_points - 1) * 32 floats
__global__
void transmissionSpectrumKernel(
  const int nb_spectral_points,
  const int nb_grid_points,
  const float bottom_radius,
  const float star_radius,
  const float* __restrict__ altitude_ns,
  const float* __restrict__ path_lengths,
  const float* __restrict__ absorption_coeff,
  const float* __restrict__ scattering_coeff,
  const float* __restrict__ cloud_optical_depth,
  float* spectrum)
{ 
  extern __shared__ float shared[];

  //extinction coefficient per level, later re-used for the 
  //transmission[t][wavelength of the block] of the tangent heights t
  float* transmission = shared;
  //sum of the extinction coefficients of the two levels enclosing a layer
  float* layer_extinction = shared + nb_grid_points * transmission_block_wavelengths;

  const int w = blockIdx.x * transmission_block_wavelengths + threadIdx.x;
  const bool active = w < nb_spectral_points;

  for (int i = threadIdx.y; i < nb_grid_points; i += transmission_block_heights)
  {
    float extinction = 0.0f;

    if (active)
    {
      const size_t j = static_cast<size_t>(i)*nb_spectral_points + w;

      extinction = absorption_coeff[j];

      if (scattering_coeff != nullptr)
        extinction += scattering_coeff[j];

      if (cloud_optical_depth != nullptr && i > 0)
        extinction += cloud_optical_depth[j - nb_spectral_points]
                    / (altitude_ns[i] - altitude_ns[i-1]);
    }

    transmission[i * transmission_block_wavelengths + threadIdx.x] = extinction;
  }

  __syncthreads();

  for (int i = threadIdx.y; i < nb_grid_points-1; i += transmission_block_heights)
    layer_extinction[i * transmission_block_wavelengths + threadIdx.x] =
        transmission[i * transmission_block_wavelengths + threadIdx.x]
      + transmission[(i+1) * transmission_block_wavelengths + threadIdx.x];

  __syncthreads();

  for (int t = nb_grid_points - 2 - threadIdx.y; t >= 0; t -= transmission_block_heights)
  {
    float tau = 0.0f;

    if (active)
    {
      const float* path = path_lengths + t * nb_grid_points;

      for (int i = t; i < nb_grid_points-1; ++i)
      {
        tau += path[i] * layer_extinction[i * transmission_block_wavelengths + threadIdx.x];

        if (tau > static_cast<float>(transmission_optical_depth_cutoff)) break;
      }
    }

    transmission[t * transmission_block_wavelengths + threadIdx.x] = __expf(-tau);
  }

  __syncthreads();

  if (threadIdx.y != 0 || !active) return;

  //integrate effective tangent height
  //single precision throughout: double-precision operations are up to 64 times
  //slower on consumer GPUs and would dominate the run time of the kernel
  float effective_tangent_height = 0.0f;

  float prev_transmission = 1.0f;  // top limb = transparent
  float prev_alt = altitude_ns[nb_grid_points-1];

  for (int t = nb_grid_points-2; t >= 0; --t)
  {
    const float transmission_t = transmission[t * transmission_block_wavelengths + threadIdx.x];

    //trapezoidal integration for effective tangent height
    //the factor of 0.5 from the trapezoidal rule cancels with 
    //the factor of 2 from accounting for both hemispheres
    const float alt = altitude_ns[t];

    effective_tangent_height +=
      ( (bottom_radius + alt) * (1.0f - transmission_t)
      + (bottom_radius + prev_alt) * (1.0f - prev_transmission) )
      * (prev_alt - alt);

    prev_transmission = transmission_t;
    prev_alt = alt;

    //if the optical depth is too high, we simply add the bottom part
    if (transmission_t < static_cast<float>(transmission_cutoff))
    {
      effective_tangent_height += (2.0f*bottom_radius + alt)*alt;
      break;
    }
  }

  //transit depth in ppm: (planet radius)^2 = bottom_radius^2 + effective_tangent_height
  spectrum[w] = (effective_tangent_height + bottom_radius*bottom_radius) /
                (star_radius * star_radius) * 1e6f;
}



__host__ 
void  TransmissionModel::calcTransitDepthGPU(
  float* transit_radius_dev,
  float* absorption_coeff_dev, 
  float* scattering_coeff_dev, 
  float* cloud_optical_depth_dev, 
  const Atmosphere& atmosphere, 
  const size_t nb_spectral_points, 
  const double radius_planet, 
  const double radius_star)
{
  const int threads = 256;

  //path lengths of the limb rays for the current atmospheric structure
  const size_t path_size = nb_grid_points * nb_grid_points;

  if (transmission_path_dev == nullptr || transmission_path_size != path_size)
  {
    if (transmission_path_dev != nullptr) deleteFromDevice(transmission_path_dev);

    allocateOnDevice(transmission_path_dev, path_size);
    transmission_path_size = path_size;
  }

  const int blocks_path = (path_size + threads - 1) / threads;

  transmissionPathLengthKernel<<<blocks_path, threads>>>(
    nb_grid_points,
    static_cast<float>(radius_planet),
    atmosphere.altitude_dev,
    transmission_path_dev);

  CUDA_CHECK_AFTER_KERNEL();

  const dim3 block_dim(transmission_block_wavelengths, transmission_block_heights);
  const int blocks = (nb_spectral_points + transmission_block_wavelengths - 1) / transmission_block_wavelengths;
  const size_t shared_bytes = (2*nb_grid_points - 1) * transmission_block_wavelengths * sizeof(float);

  static size_t attribute_bytes = 0;

  if (shared_bytes > 48 * 1024 && shared_bytes > attribute_bytes)
  {
    gpuErrchk(cudaFuncSetAttribute(
      transmissionSpectrumKernel,
      cudaFuncAttributeMaxDynamicSharedMemorySize,
      shared_bytes));
    attribute_bytes = shared_bytes;
  }

  transmissionSpectrumKernel<<<blocks, block_dim, shared_bytes>>>(
    nb_spectral_points,
    nb_grid_points,
    static_cast<float>(radius_planet),
    static_cast<float>(radius_star),
    atmosphere.altitude_dev,
    transmission_path_dev,
    absorption_coeff_dev,
    scattering_coeff_dev,
    cloud_optical_depth_dev,
    transit_radius_dev);

  CUDA_CHECK_AFTER_KERNEL();
}

}