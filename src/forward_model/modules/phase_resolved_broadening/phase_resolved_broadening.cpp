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
 * quadratically limb-darkened stellar disc at mid-transit.  Each ring point in
 * front of the disc contributes at its line-of-sight velocity, weighted by the
 * stellar intensity at that position; the resulting velocity distribution is
 * convolved with the instrumental Gaussian.  The kernel is tabulated in velocity
 * and applied as a normalised convolution on the actual velocity steps of the
 * grid (see ../velocity_space_convolution.h).
 */


#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <omp.h>

#include "phase_resolved_broadening.h"

#include "../velocity_space_convolution.h"
#include "../../../spectral_grid/spectral_grid.h"
#include "../../../additional/physical_const.h"
#include "phase_resolved_convolution.h"


namespace bear{


PhaseResolvedBroadening::PhaseResolvedBroadening(
  const std::vector<std::string>& parameters,
  SpectralGrid* spectral_grid_)
  : spectral_grid(spectral_grid_)
{
  // 7 parameters for the Brogi et al. (2016) transit geometry / broadening.
  //parameter_names is the source of truth and must match the read order in
  //modifySpectrum / modifySpectrumGPU (parameter[0..6]):
  //  [0] v_eq       -> equatorial rotation velocity (km/s)
  //  [1] v_wind     -> equatorial super-rotation wind velocity (km/s)
  //  [2] sigma_inst -> instrumental Gaussian width (km/s)
  //  [3] u1         -> linear quadratic-limb-darkening coefficient
  //  [4] u2         -> quadratic quadratic-limb-darkening coefficient
  //  [5] Rp_Rs      -> planet-to-star radius ratio
  //  [6] impact_b   -> transit impact parameter (stellar radii)
  parameter_names = {
    "v_eq", "v_wind", "sigma_inst", "u1", "u2", "rp_rs", "impact_b"};

  setVelocitySteps();

  std::cout << "Phase-resolved broadening module initialised, velocity step "
            << dv_min_kms << " km/s (minimum)\n";
}


void PhaseResolvedBroadening::setVelocitySteps()
{
  dv_kms = gridVelocitySteps(spectral_grid->wavelength_list);

  dv_min_kms = dv_kms.empty() ? 0.0 : *std::min_element(dv_kms.begin(), dv_kms.end());
}


void PhaseResolvedBroadening::freeDeviceBuffers()
{
  if (temp_buffer_gpu != nullptr)
    deleteFromDevice(temp_buffer_gpu);

  if (table_gpu != nullptr)
    deleteFromDevice(table_gpu);

  table_gpu_capacity = 0;

  if (dv_gpu != nullptr)
    deleteFromDevice(dv_gpu);
}


void PhaseResolvedBroadening::setSpectralGrid(SpectralGrid* grid)
{
  spectral_grid = grid;

  freeDeviceBuffers();

  setVelocitySteps();

  const double dv_max = dv_kms.empty() ? 0.0 : *std::max_element(dv_kms.begin(), dv_kms.end());

  std::cout << "Phase-resolved broadening: updated spectral grid, velocity step "
            << dv_min_kms << " - " << dv_max << " km/s, "
            << spectral_grid->wavelength_list.size() << " spectral points\n";
}


PhaseResolvedBroadening::~PhaseResolvedBroadening()
{
  freeDeviceBuffers();
}


void PhaseResolvedBroadening::buildBroadeningKernel(
  double v_eq, double v_wind, double sigma_inst,
  double u1, double u2, double Rp_Rs, double impact_b,
  std::vector<double>& table, double& v_min, double& step, double& half_width)
{
  // A zero-width instrumental profile would make the Gaussian singular
  const double sigma = std::max(std::fabs(sigma_inst), 1e-3);

  // Maximum possible velocity extent + 5-sigma Gaussian wings
  half_width = std::fabs(v_eq) + std::fabs(v_wind) + 5.0 * sigma;

  // Fine enough to resolve both the kernel (width ~sigma) and the grid steps it is
  // sampled at, so that interpolating the table costs well below 1e-3 in accuracy.
  step = std::max(std::min(sigma / 8.0, dv_min_kms / 4.0), 1e-4);

  const int n_half = static_cast<int>(std::ceil(half_width / step));
  const int size = 2 * n_half + 1;
  v_min = -n_half * step;

  // Distribution of the ring's line-of-sight velocities, weighted by the local
  // stellar intensity (cloud-in-cell assignment to the table points)
  std::vector<double> histogram(size, 0.0);

  const double sin_wind_limit = std::sin(WIND_LAT_LIMIT_DEG * constants::pi / 180.0);
  const double dtheta = 2.0 * constants::pi / N_ANGLE;

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

    const double x = (v_total - v_min) / step;
    const int k = static_cast<int>(std::floor(x));
    const double t = x - k;

    if (k >= 0 && k < size)
      histogram[k] += I_mu * (1.0 - t);

    if (k + 1 >= 0 && k + 1 < size)
      histogram[k + 1] += I_mu * t;
  }

  // Convolve the velocity distribution with the instrumental Gaussian.  The
  // overall normalisation is irrelevant: the convolution normalises per point.
  const int g_half = static_cast<int>(std::ceil(5.0 * sigma / step));
  std::vector<double> gauss(2 * g_half + 1);

  for (int g = -g_half; g <= g_half; ++g)
  {
    const double v = g * step;
    gauss[g + g_half] = std::exp(-v * v / (2.0 * sigma * sigma));
  }

  table.assign(size, 0.0);

  for (int k = 0; k < size; ++k)
  {
    if (histogram[k] == 0.0)
      continue;

    for (int g = -g_half; g <= g_half; ++g)
    {
      const int j = k + g;

      if (j >= 0 && j < size)
        table[j] += histogram[k] * gauss[g + g_half];
    }
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

  std::vector<double> table;
  double v_min, step, half_width;
  buildBroadeningKernel(v_eq, v_wind, sigma_inst, u1, u2, Rp_Rs, impact_b,
                        table, v_min, step, half_width);

  // A kernel narrower than every grid step leaves the spectrum unchanged
  if (half_width < dv_min_kms)
    return;

  std::vector<double> broadened;
  convolveVelocitySpace(spectrum, broadened, dv_kms, half_width,
                        TabulatedProfile(table, v_min, step));

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

  // Build the kernel table on the CPU
  std::vector<double> table;
  double v_min, step, half_width;
  buildBroadeningKernel(v_eq, v_wind, sigma_inst, u1, u2, Rp_Rs, impact_b,
                        table, v_min, step, half_width);

  if (half_width < dv_min_kms)
    return;

  const int n_pixels = static_cast<int>(spectral_grid->nbSpectralPoints());

  if (temp_buffer_gpu == nullptr)
    allocateOnDevice(temp_buffer_gpu, static_cast<size_t>(n_pixels));

  if (dv_gpu == nullptr)
  {
    std::vector<float> dv_f(dv_kms.begin(), dv_kms.end());
    moveToDevice(dv_gpu, dv_f);
  }

  // Upload the table, growing the buffer with some headroom when needed
  if (table_gpu == nullptr || table_gpu_capacity < table.size())
  {
    if (table_gpu != nullptr)
      deleteFromDevice(table_gpu);

    table_gpu_capacity = std::max(table.size(), static_cast<size_t>(1024));
    allocateOnDevice(table_gpu, table_gpu_capacity);
  }

  std::vector<float> table_f(table.begin(), table.end());
  moveToDevice(table_gpu, table_f);

  // Copy spectrum to temp buffer, convolve temp -> spectrum_gpu
  copyOnDevice(temp_buffer_gpu, spectrum_gpu, static_cast<size_t>(n_pixels));

  applyTabulatedConvolutionGPU(
    temp_buffer_gpu, spectrum_gpu, n_pixels, dv_gpu,
    table_gpu, static_cast<int>(table.size()), v_min, step, half_width);
}


}
