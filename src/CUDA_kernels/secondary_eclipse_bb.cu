// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <iostream>
#include <cmath>
#include <stdio.h>

#include "../forward_model/secondary_eclipse_bb/secondary_eclipse_bb.h"

#include "error_check.h"
#include "planck_function.h"
#include "../additional/physical_const.h"


namespace bear{


__global__ void planetBlackBodyFlux(
  double* wavenumbers,
  int nb_wavenumbers,
  const double planet_temperature,
  float* flux)
{
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_wavenumbers; tid += blockDim.x * gridDim.x)
  {
    flux[tid] = static_cast<float>(planckFunction(
      planet_temperature,
      wavenumbers[tid]) * constants::pi * 1e-3);
  }
}



__global__ void OccultationBBDevice(
  float* secondary_eclipse,
  float* planet_spectrum,
  const float* stellar_spectrum,
  const int nb_points,
  const double radius_ratio,
  const float* albedo_contribution)
{

  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < nb_points; i += blockDim.x * gridDim.x)
  {
    secondary_eclipse[i] = static_cast<float>(static_cast<double>(planet_spectrum[i])/static_cast<double>(stellar_spectrum[i]) * radius_ratio*radius_ratio * 1e6);

    if (albedo_contribution != nullptr)
      secondary_eclipse[i] += static_cast<float>(static_cast<double>(albedo_contribution[i])*1e6);
  }
}


__host__ void OccultationBlackBodyModel::calcPlanetSpectrumGPU(
  const double planet_temperature,
  float* spectrum_dev)
{
  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  planetBlackBodyFlux<<<blocks,threads>>>(
    spectral_grid->wavenumber_list_gpu,
    nb_spectral_points,
    planet_temperature,
    spectrum_dev);

  CUDA_CHECK_AFTER_KERNEL();
}


__host__ void OccultationBlackBodyModel::calcOccultationGPU(
  float* secondary_eclipse,
  float* planet_spectrum,
  const float* stellar_spectrum,
  const int nb_points,
  const double radius_ratio,
  const float* albedo_contribution)
{
  int threads = 256;

  int blocks = nb_points / threads;
  if (nb_points % threads) blocks++;


  OccultationBBDevice<<<blocks,threads>>>(
    secondary_eclipse,
    planet_spectrum,
    stellar_spectrum,
    nb_points,
    radius_ratio,
    albedo_contribution);


  CUDA_CHECK_AFTER_KERNEL();
}


}