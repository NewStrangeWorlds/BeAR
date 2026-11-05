// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <iostream>
#include <cmath>
#include <stdio.h>

#include "../forward_model/stellar_spectrum/stellar_spectrum_grid.h"

#include "error_check.h"
#include "planck_function.h"
#include "../additional/physical_const.h"


namespace bear{


__global__
void stellarSpectrumInterpolationFl(
  const float xd, const float yd, const float zd,
  const float* __restrict__ c000, const float* __restrict__ c100,
  const float* __restrict__ c010, const float* __restrict__ c110,
  const float* __restrict__ c001, const float* __restrict__ c101,
  const float* __restrict__ c011, const float* __restrict__ c111,
  const int nb_wavenumbers,
  float* __restrict__ spectrum_dev)
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
  float* spectrum_gpu)
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

  CUDA_CHECK_AFTER_KERNEL();
}


}
 
