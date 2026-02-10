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
#include <cmath>
#include <stdio.h>

#include "../forward_model/stellar_spectrum/stellar_spectrum_grid.h"

#include "error_check.h"
#include "planck_function.h"
#include "../additional/physical_const.h"


namespace bear{


__global__ void stellarSpectrumInterpolationOld(
  const double xd, const double yd, const double zd,
  const double* c000, const double* c100,
  const double* c010, const double* c110,
  const double* c001, const double* c101,
  const double* c011, const double* c111,
  const int nb_wavenumbers,
  double* spectrum_dev)
{
  for (int tid = blockIdx.x * blockDim.x + threadIdx.x; tid < nb_wavenumbers; tid += blockDim.x * gridDim.x)
  {
    const double c00 = c000[tid] * (1. - xd) + c100[tid]*xd;
    const double c01 = c001[tid] * (1. - xd) + c101[tid]*xd;
    const double c10 = c010[tid] * (1. - xd) + c110[tid]*xd;
    const double c11 = c011[tid] * (1. - xd) + c111[tid]*xd;

    const double c0 = c00 * (1. - yd) + c10*yd;
    const double c1 = c01 * (1. - yd) + c11*yd;

    const double c = c0 * (1. - zd) + c1 * zd;

    spectrum_dev[tid] = c0 * (1. - zd) + c1 * zd;
  }
}


__global__ 
void stellarSpectrumInterpolation(
  const double xd, const double yd, const double zd,
  const double* __restrict__ c000, const double* __restrict__ c100,
  const double* __restrict__ c010, const double* __restrict__ c110,
  const double* __restrict__ c001, const double* __restrict__ c101,
  const double* __restrict__ c011, const double* __restrict__ c111,
  const int nb_wavenumbers,
  double* __restrict__ spectrum_dev)
{
  // Pre-calculate weights to save subtractions in the loop
  const double omx = 1.0 - xd;
  const double omy = 1.0 - yd;
  const double omz = 1.0 - zd;

  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;

  for (; tid < nb_wavenumbers; tid += stride)
  {
    // Layer 0 (z=0)
    double xy0 = omy * (c000[tid] * omx + c100[tid] * xd) + 
                 yd  * (c010[tid] * omx + c110[tid] * xd);
    
    // Layer 1 (z=1)
    double xy1 = omy * (c001[tid] * omx + c101[tid] * xd) + 
                 yd  * (c011[tid] * omx + c111[tid] * xd);

    // Final Interpolation
    spectrum_dev[tid] = xy0 * omz + xy1 * zd;
  }
}


__global__ 
void stellarSpectrumInterpolationFl(
  const float xd, const float yd, const float zd,
  const float* __restrict__ c000, const float* __restrict__ c100,
  const float* __restrict__ c010, const float* __restrict__ c110,
  const float* __restrict__ c001, const float* __restrict__ c101,
  const float* __restrict__ c011, const float* __restrict__ c111,
  const int nb_wavenumbers,
  double* __restrict__ spectrum_dev)
{
  // Pre-calculate weights to save subtractions in the loop
  const float omx = 1.0f - xd;
  const float omy = 1.0f - yd;
  const float omz = 1.0f - zd;

  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;

  for (; tid < nb_wavenumbers; tid += stride)
  {
    // Layer 0 (z=0)
    float xy0 = omy * (c000[tid] * omx + c100[tid] * xd) + 
                 yd  * (c010[tid] * omx + c110[tid] * xd);
    
    // Layer 1 (z=1)
    float xy1 = omy * (c001[tid] * omx + c101[tid] * xd) + 
                 yd  * (c011[tid] * omx + c111[tid] * xd);

    // Final Interpolation
    spectrum_dev[tid] = xy0 * omz + xy1 * zd;
  }
}



__host__ void StellarSpectrumGrid::calcFluxGPU(
  const std::vector<double>& parameter,
  double* spectrum_gpu)
{
  double temperature_int = parameter[0];
  double log_g_int = parameter[1];
  double metallicity_int = parameter[2];

  if (temperature_int < effective_temperature.front()) temperature_int = effective_temperature.front();
  if (temperature_int > effective_temperature.back()) temperature_int = effective_temperature.back();

  if (log_g_int < log_g.front()) log_g_int = log_g.front();
  if (log_g_int > log_g.back()) log_g_int = log_g.back();

  if (metallicity_int < metallicity.front()) metallicity_int = metallicity.front();
  if (metallicity_int > metallicity.back()) metallicity_int = metallicity.back();


  std::pair<size_t, size_t> t_i = findClosestGridPoints(effective_temperature, temperature_int);
  std::pair<size_t, size_t> g_i = findClosestGridPoints(log_g, log_g_int);
  std::pair<size_t, size_t> m_i = findClosestGridPoints(metallicity, metallicity_int);

  grid[t_i.first][g_i.first][m_i.first]->sampleSpectrum(spectral_grid, grid_wavelengths, true);
  grid[t_i.second][g_i.first][m_i.first]->sampleSpectrum(spectral_grid, grid_wavelengths, true);
  grid[t_i.first][g_i.second][m_i.first]->sampleSpectrum(spectral_grid, grid_wavelengths, true);
  grid[t_i.second][g_i.second][m_i.first]->sampleSpectrum(spectral_grid, grid_wavelengths, true);
  grid[t_i.first][g_i.first][m_i.second]->sampleSpectrum(spectral_grid, grid_wavelengths, true);
  grid[t_i.second][g_i.first][m_i.second]->sampleSpectrum(spectral_grid, grid_wavelengths, true);
  grid[t_i.first][g_i.second][m_i.second]->sampleSpectrum(spectral_grid, grid_wavelengths, true);
  grid[t_i.second][g_i.second][m_i.second]->sampleSpectrum(spectral_grid, grid_wavelengths, true);

  double xd = (temperature_int - effective_temperature[t_i.first])/(effective_temperature[t_i.second] - effective_temperature[t_i.first]);
  if (temperature_int == effective_temperature[t_i.first]) xd = 0;

  double yd = (log_g_int - log_g[g_i.first])/(log_g[g_i.second] - log_g[g_i.first]);
  if (log_g_int == log_g[g_i.first]) yd = 0;

  double zd = (metallicity_int - metallicity[m_i.first])/(metallicity[m_i.second] - metallicity[m_i.first]);
  if (metallicity_int == metallicity[m_i.first]) zd = 0;


  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();
  int threads = 256;
  int blocks = nb_spectral_points / threads;
  if (nb_spectral_points % threads) blocks++;

  const double effective_temperature = parameter[0];

  stellarSpectrumInterpolationFl<<<blocks,threads>>>(
    xd, yd, zd,
    grid[t_i.first][g_i.first][m_i.first]->spectrum_gpu,    //c000
    grid[t_i.second][g_i.first][m_i.first]->spectrum_gpu,   //c100
    grid[t_i.first][g_i.second][m_i.first]->spectrum_gpu,   //c010
    grid[t_i.second][g_i.second][m_i.first]->spectrum_gpu,  //c110
    grid[t_i.first][g_i.first][m_i.second]->spectrum_gpu,   //c001
    grid[t_i.second][g_i.first][m_i.second]->spectrum_gpu,  //c101
    grid[t_i.first][g_i.second][m_i.second]->spectrum_gpu,  //c011
    grid[t_i.second][g_i.second][m_i.second]->spectrum_gpu, //c111
    nb_spectral_points,
    spectrum_gpu);

  cudaDeviceSynchronize();
  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );
}


}
 
