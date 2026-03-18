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
 * Phase-resolved broadening module based on Brogi et al. (2016) Section 4.2.
 *
 * Constructs a 2D model of the planet's atmospheric ring transiting a
 * quadratically limb-darkened stellar disc at mid-transit.  Each ring pixel
 * contributes a Gaussian (instrumental profile) centered on its line-of-sight
 * velocity, weighted by the stellar intensity at that position.  The sum
 * produces a 1D broadening kernel that is convolved with the spectrum.
 */


#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <omp.h>

#include "phase_resolved_broadening.h"

#include "../../../spectral_grid/spectral_grid.h"
#include "../../../additional/physical_const.h"
#include "phase_resolved_convolution.h"


namespace bear{


PhaseResolvedBroadening::PhaseResolvedBroadening(
  const std::vector<std::string>& parameters,
  SpectralGrid* spectral_grid_)
  : spectral_grid(spectral_grid_)
{
  nb_parameters = 7;

  const size_t n = spectral_grid->wavelength_list.size();

  if (n >= 2)
  {
    const size_t mid = n / 2;
    const double dlambda = std::fabs(spectral_grid->wavelength_list[mid+1]
                                    - spectral_grid->wavelength_list[mid]);
    const double lambda = spectral_grid->wavelength_list[mid];
    delta_v_kms = (constants::light_c * 1e-5) * dlambda / lambda;
  }

  std::cout << "Phase-resolved broadening module initialised, delta_v = "
            << delta_v_kms << " km/s\n";
}


void PhaseResolvedBroadening::setSpectralGrid(SpectralGrid* grid)
{
  spectral_grid = grid;

  if (temp_buffer_gpu != nullptr)
  {
    deleteFromDevice(temp_buffer_gpu);
    temp_buffer_gpu = nullptr;
  }

  if (kernel_gpu != nullptr)
  {
    deleteFromDevice(kernel_gpu);
    kernel_gpu = nullptr;
    kernel_gpu_size = 0;
  }

  const size_t n = spectral_grid->wavelength_list.size();

  if (n >= 2)
  {
    const size_t mid = n / 2;
    const double dlambda = std::fabs(spectral_grid->wavelength_list[mid+1]
                                    - spectral_grid->wavelength_list[mid]);
    const double lambda = spectral_grid->wavelength_list[mid];
    delta_v_kms = (constants::light_c * 1e-5) * dlambda / lambda;
  }

  std::cout << "Phase-resolved broadening: updated spectral grid, delta_v = "
            << delta_v_kms << " km/s, " << n << " spectral points\n";
}


PhaseResolvedBroadening::~PhaseResolvedBroadening()
{
  if (temp_buffer_gpu != nullptr)
    deleteFromDevice(temp_buffer_gpu);

  if (kernel_gpu != nullptr)
    deleteFromDevice(kernel_gpu);
}


void PhaseResolvedBroadening::buildBroadeningKernel(
  double v_eq, double v_wind, double sigma_inst,
  double u1, double u2, double Rp_Rs, double impact_b,
  std::vector<double>& kernel_out, int& kernel_hw)
{
  // Maximum possible velocity extent + 5-sigma Gaussian wings
  const double v_max = std::fabs(v_eq) + std::fabs(v_wind)
                     + 5.0 * sigma_inst;
  kernel_hw = static_cast<int>(std::ceil(v_max / delta_v_kms));

  if (kernel_hw < 1)
    kernel_hw = 1;

  const int kernel_size = 2 * kernel_hw + 1;
  kernel_out.assign(kernel_size, 0.0);

  const double inv_2sig2 = 1.0 / (2.0 * sigma_inst * sigma_inst);
  const double gauss_norm = 1.0 / (sigma_inst * std::sqrt(2.0 * constants::pi));
  const double sin_wind_limit = std::sin(WIND_LAT_LIMIT_DEG * constants::pi / 180.0);
  const double dtheta = 2.0 * constants::pi / N_ANGLE;

  double total_weight = 0.0;

  for (int a = 0; a < N_ANGLE; ++a)
  {
    const double theta = a * dtheta;
    const double cos_theta = std::cos(theta);
    const double sin_theta = std::sin(theta);

    // Ring pixel position in stellar radii
    const double x_ring = Rp_Rs * cos_theta;
    const double y_ring = impact_b + Rp_Rs * sin_theta;

    // Check if on stellar disc
    const double r2 = x_ring * x_ring + y_ring * y_ring;
    if (r2 >= 1.0)
      continue;

    // Quadratic limb darkening: I(mu) = 1 - u1*(1-mu) - u2*(1-mu)^2
    const double mu = std::sqrt(1.0 - r2);
    const double one_minus_mu = 1.0 - mu;
    const double I_mu = 1.0 - u1 * one_minus_mu - u2 * one_minus_mu * one_minus_mu;

    if (I_mu <= 0.0)
      continue;

    // Line-of-sight velocity from rigid-body rotation
    double v_total = v_eq * cos_theta;

    // Equatorial super-rotation wind (within +/-25 deg latitude)
    if (std::fabs(sin_theta) < sin_wind_limit)
    {
      // Receding limb (x > 0): wind adds velocity
      // Approaching limb (x < 0): wind subtracts velocity
      if (x_ring > 0.0)
        v_total += v_wind;
      else if (x_ring < 0.0)
        v_total -= v_wind;
    }

    // Accumulate this ring pixel's Gaussian contribution to the kernel
    for (int k = 0; k < kernel_size; ++k)
    {
      const double v_pixel = (k - kernel_hw) * delta_v_kms;
      const double dv = v_pixel - v_total;
      const double g = gauss_norm * std::exp(-dv * dv * inv_2sig2);
      kernel_out[k] += I_mu * g;
    }

    total_weight += I_mu;
  }

  // Normalize kernel to unit sum (flux conservation)
  double kernel_sum = 0.0;
  for (int k = 0; k < kernel_size; ++k)
    kernel_sum += kernel_out[k];

  if (kernel_sum > 0.0)
  {
    const double inv_sum = 1.0 / kernel_sum;
    for (int k = 0; k < kernel_size; ++k)
      kernel_out[k] *= inv_sum;
  }
}


void PhaseResolvedBroadening::convolveWithKernelCPU(
  const std::vector<double>& spectrum_in,
  std::vector<double>& spectrum_out,
  const std::vector<double>& kernel,
  int kernel_hw)
{
  const int n = static_cast<int>(spectrum_in.size());
  spectrum_out.resize(n);

  #pragma omp parallel for schedule(dynamic, 64)
  for (int i = 0; i < n; ++i)
  {
    double sum = 0.0;
    const int j0 = std::max(0, i - kernel_hw);
    const int j1 = std::min(n - 1, i + kernel_hw);

    for (int j = j0; j <= j1; ++j)
      sum += kernel[j - i + kernel_hw] * spectrum_in[j];

    spectrum_out[i] = sum;
  }
}


void PhaseResolvedBroadening::modifySpectrum(
  const std::vector<double>& parameter,
  Atmosphere* atmosphere,
  std::vector<double>& spectrum)
{
  const double v_eq       = parameter[0];
  const double v_wind     = parameter[1];
  const double sigma_inst = parameter[2];
  const double u1         = parameter[3];
  const double u2         = parameter[4];
  const double Rp_Rs      = parameter[5];
  const double impact_b   = parameter[6];

  // Nothing to broaden if all velocities are sub-pixel
  if (std::fabs(v_eq) + std::fabs(v_wind) < 0.01 * delta_v_kms
      && sigma_inst < 0.01 * delta_v_kms)
    return;

  std::vector<double> kernel;
  int kernel_hw;
  buildBroadeningKernel(v_eq, v_wind, sigma_inst, u1, u2, Rp_Rs, impact_b,
                        kernel, kernel_hw);

  std::vector<double> broadened;
  convolveWithKernelCPU(spectrum, broadened, kernel, kernel_hw);

  spectrum = std::move(broadened);
}


void PhaseResolvedBroadening::modifySpectrumGPU(
  const std::vector<double>& parameter,
  Atmosphere* atmosphere,
  float* spectrum_gpu)
{
  const double v_eq       = parameter[0];
  const double v_wind     = parameter[1];
  const double sigma_inst = parameter[2];
  const double u1         = parameter[3];
  const double u2         = parameter[4];
  const double Rp_Rs      = parameter[5];
  const double impact_b   = parameter[6];

  if (std::fabs(v_eq) + std::fabs(v_wind) < 0.01 * delta_v_kms
      && sigma_inst < 0.01 * delta_v_kms)
    return;

  // Build kernel on CPU
  std::vector<double> kernel;
  int kernel_hw;
  buildBroadeningKernel(v_eq, v_wind, sigma_inst, u1, u2, Rp_Rs, impact_b,
                        kernel, kernel_hw);

  const int n_pixels = static_cast<int>(spectral_grid->nbSpectralPoints());
  const int kernel_size = 2 * kernel_hw + 1;

  // Allocate GPU buffers
  if (temp_buffer_gpu == nullptr)
    allocateOnDevice(temp_buffer_gpu, static_cast<size_t>(n_pixels));

  if (kernel_gpu == nullptr || kernel_gpu_size < kernel_size)
  {
    if (kernel_gpu != nullptr)
      deleteFromDevice(kernel_gpu);

    // Allocate with some headroom to avoid frequent reallocation
    kernel_gpu_size = std::max(kernel_size, 512);
    allocateOnDevice(kernel_gpu, static_cast<size_t>(kernel_gpu_size));
  }

  // Upload kernel to GPU (convert double -> float)
  std::vector<float> kernel_f(kernel_size);
  for (int k = 0; k < kernel_size; ++k)
    kernel_f[k] = static_cast<float>(kernel[k]);

  moveToDevice(kernel_gpu, kernel_f);

  // Copy spectrum to temp buffer, convolve temp -> spectrum_gpu
  copyOnDevice(temp_buffer_gpu, spectrum_gpu, static_cast<size_t>(n_pixels));

  applyPrecomputedConvolutionGPU(
    temp_buffer_gpu, spectrum_gpu, kernel_gpu, kernel_hw, n_pixels);
}


}
