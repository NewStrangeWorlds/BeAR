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


#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <iomanip>

#include "velocity_broadening.h"

#include "../../../spectral_grid/spectral_grid.h"
#include "../../../additional/physical_const.h"
#include "../../../additional/aux_functions.h"
#include "../../../additional/exceptions.h"
#include "../../../CUDA_kernels/highres_convolution.h"


namespace bear{


VelocityBroadening::VelocityBroadening (
  const std::vector<std::string>& velocity_broadening_parameters,
  SpectralGrid* spectral_grid_)
  : spectral_grid(spectral_grid_)
{
  // 4 parameters: vsini (km/s), sigma_inst (km/s), epsilon, v_wind (km/s)
  //parameter_names is the source of truth and must match the read order in
  //modifySpectrum (parameter[0..3]).
  parameter_names = {"vsini", "sigma_inst", "epsilon", "v_wind"};

  // Compute velocity spacing from spectral grid.
  // On a log-lambda grid: delta_v = c * (lambda[i+1] - lambda[i]) / lambda[i]
  const size_t n = spectral_grid->wavelength_list.size();

  if (n >= 2)
  {
    const size_t mid = n / 2;
    const double dlambda = std::fabs(spectral_grid->wavelength_list[mid+1]
                                    - spectral_grid->wavelength_list[mid]);
    const double lambda = spectral_grid->wavelength_list[mid];
    delta_v_kms = (constants::light_c * 1e-5) * dlambda / lambda;
  }

  std::cout << "Velocity broadening module initialised, delta_v = "
            << delta_v_kms << " km/s\n";
}


void VelocityBroadening::setSpectralGrid(SpectralGrid* grid)
{
  spectral_grid = grid;

  // Free old GPU buffers since size may change
  if (temp_buffer_gpu != nullptr)
  {
    deleteFromDevice(temp_buffer_gpu);
    temp_buffer_gpu = nullptr;
  }

  if (temp_buffer2_gpu != nullptr)
  {
    deleteFromDevice(temp_buffer2_gpu);
    temp_buffer2_gpu = nullptr;
  }

  // Recompute velocity spacing from new grid
  const size_t n = spectral_grid->wavelength_list.size();

  if (n >= 2)
  {
    const size_t mid = n / 2;
    const double dlambda = std::fabs(spectral_grid->wavelength_list[mid+1]
                                    - spectral_grid->wavelength_list[mid]);
    const double lambda = spectral_grid->wavelength_list[mid];
    delta_v_kms = (constants::light_c * 1e-5) * dlambda / lambda;
  }

  std::cout << "Velocity broadening: updated spectral grid, delta_v = "
            << delta_v_kms << " km/s, " << n << " spectral points\n";
}


VelocityBroadening::~VelocityBroadening()
{
  if (temp_buffer_gpu != nullptr)
    deleteFromDevice(temp_buffer_gpu);

  if (temp_buffer2_gpu != nullptr)
    deleteFromDevice(temp_buffer2_gpu);
}


static double rotKernel(double dx, double vsini_pix, double epsilon)
{
  double x = dx / vsini_pix;
  if (std::fabs(x) >= 1.0)
    return 0.0;

  double one_minus_x2 = 1.0 - x * x;
  double c1 = 2.0 * (1.0 - epsilon);
  double c2 = 0.5 * constants::pi * epsilon;
  double denom = constants::pi * vsini_pix * (1.0 - epsilon / 3.0);

  return (c1 * std::sqrt(one_minus_x2) + c2 * one_minus_x2) / denom;
}


static double gaussKernel(double dx, double sigma_pix)
{
  double inv_2sig2 = 1.0 / (2.0 * sigma_pix * sigma_pix);
  double norm = 1.0 / (sigma_pix * std::sqrt(2.0 * constants::pi));
  return norm * std::exp(-dx * dx * inv_2sig2);
}


void VelocityBroadening::convolveSpectrumCPU(
  const std::vector<double>& spectrum_in,
  std::vector<double>& spectrum_out,
  double sigma_kms,
  double vsini_kms,
  double epsilon)
{
  const int n = static_cast<int>(spectrum_in.size());
  spectrum_out.resize(n);

  const double sigma_pix = sigma_kms / delta_v_kms;
  const int gauss_hw = static_cast<int>(std::ceil(5.0 * sigma_pix));

  const bool do_rotation = vsini_kms > 0.5 * delta_v_kms;
  const bool do_gaussian = sigma_pix > 0.01;

  if (do_rotation && do_gaussian)
  {
    const double vsini_pix = vsini_kms / delta_v_kms;
    const int rot_hw = static_cast<int>(std::ceil(vsini_pix));

    std::vector<double> temp(n);

    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < n; ++i)
    {
      double sum = 0.0;
      const int j0 = std::max(0, i - rot_hw);
      const int j1 = std::min(n - 1, i + rot_hw);

      for (int j = j0; j <= j1; ++j)
        sum += rotKernel(static_cast<double>(j - i), vsini_pix, epsilon)
             * spectrum_in[j];

      temp[i] = sum;
    }

    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < n; ++i)
    {
      double sum = 0.0;
      const int j0 = std::max(0, i - gauss_hw);
      const int j1 = std::min(n - 1, i + gauss_hw);

      for (int j = j0; j <= j1; ++j)
        sum += gaussKernel(static_cast<double>(j - i), sigma_pix) * temp[j];

      spectrum_out[i] = sum;
    }
  }
  else if (do_rotation)
  {
    const double vsini_pix = vsini_kms / delta_v_kms;
    const int rot_hw = static_cast<int>(std::ceil(vsini_pix));

    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < n; ++i)
    {
      double sum = 0.0;
      const int j0 = std::max(0, i - rot_hw);
      const int j1 = std::min(n - 1, i + rot_hw);

      for (int j = j0; j <= j1; ++j)
        sum += rotKernel(static_cast<double>(j - i), vsini_pix, epsilon)
             * spectrum_in[j];

      spectrum_out[i] = sum;
    }
  }
  else if (do_gaussian)
  {
    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < n; ++i)
    {
      double sum = 0.0;
      const int j0 = std::max(0, i - gauss_hw);
      const int j1 = std::min(n - 1, i + gauss_hw);

      for (int j = j0; j <= j1; ++j)
        sum += gaussKernel(static_cast<double>(j - i), sigma_pix)
             * spectrum_in[j];

      spectrum_out[i] = sum;
    }
  }
  else
  {
    spectrum_out = spectrum_in;
  }
}


void VelocityBroadening::modifySpectrum(
  const std::vector<double>& parameter,
  Atmosphere* atmosphere,
  std::vector<double>& spectrum)
{
  const double vsini_kms  = parameter[0];
  const double sigma_inst = parameter[1];
  const double epsilon    = parameter[2];
  const double v_wind     = parameter[3];

  const double sigma_total = std::sqrt(sigma_inst * sigma_inst
                                     + v_wind * v_wind);

  // Nothing to broaden — leave spectrum unchanged
  if (sigma_total < 0.01 * delta_v_kms && vsini_kms < 0.5 * delta_v_kms)
    return;

  std::vector<double> broadened;
  convolveSpectrumCPU(spectrum, broadened, sigma_total, vsini_kms, epsilon);

  spectrum = std::move(broadened);
}


void VelocityBroadening::modifySpectrumGPU(
  const std::vector<double>& parameter,
  Atmosphere* atmosphere,
  float* spectrum_gpu)
{
  const double vsini_kms  = parameter[0];
  const double sigma_inst = parameter[1];
  const double epsilon    = parameter[2];
  const double v_wind     = parameter[3];

  const double sigma_total = std::sqrt(sigma_inst * sigma_inst
                                     + v_wind * v_wind);

  // Nothing to broaden — leave spectrum unchanged
  if (sigma_total < 0.01 * delta_v_kms && vsini_kms < 0.5 * delta_v_kms)
    return;

  const int n_pixels = static_cast<int>(spectral_grid->nbSpectralPoints());

  if (temp_buffer_gpu == nullptr)
    allocateOnDevice(temp_buffer_gpu, static_cast<size_t>(n_pixels));

  if (temp_buffer2_gpu == nullptr)
    allocateOnDevice(temp_buffer2_gpu, static_cast<size_t>(n_pixels));

  // Copy input to temp buffer, then convolve temp -> spectrum_gpu (in-place)
  copyOnDevice(temp_buffer_gpu, spectrum_gpu, static_cast<size_t>(n_pixels));

  applyHighResConvolutionGPU(
    temp_buffer_gpu,
    spectrum_gpu,
    n_pixels,
    sigma_total,
    vsini_kms,
    delta_v_kms,
    epsilon,
    temp_buffer2_gpu);
}


}
