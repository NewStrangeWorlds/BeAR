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
#include <algorithm>
#include <omp.h>
#include <iomanip>

#include "velocity_broadening.h"

#include "../velocity_space_convolution.h"
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

  setVelocitySteps();

  std::cout << "Velocity broadening module initialised, velocity step "
            << dv_min_kms << " km/s (minimum)\n";
}


// Velocity steps of the current grid.  The high-res grid is a subset of the
// opacity wavenumber points, so the step varies along the grid; the minimum is
// used to skip convolutions narrower than every step (which would be exact no-ops).
void VelocityBroadening::setVelocitySteps()
{
  dv_kms = gridVelocitySteps(spectral_grid->wavelength_list);

  dv_min_kms = dv_kms.empty() ? 0.0 : *std::min_element(dv_kms.begin(), dv_kms.end());
}


void VelocityBroadening::freeDeviceBuffers()
{
  if (temp_buffer_gpu != nullptr)
    deleteFromDevice(temp_buffer_gpu);

  if (temp_buffer2_gpu != nullptr)
    deleteFromDevice(temp_buffer2_gpu);

  if (dv_gpu != nullptr)
    deleteFromDevice(dv_gpu);
}


void VelocityBroadening::setSpectralGrid(SpectralGrid* grid)
{
  spectral_grid = grid;

  // Free old GPU buffers since size may change
  freeDeviceBuffers();

  setVelocitySteps();

  const double dv_max = dv_kms.empty() ? 0.0 : *std::max_element(dv_kms.begin(), dv_kms.end());

  std::cout << "Velocity broadening: updated spectral grid, velocity step "
            << dv_min_kms << " - " << dv_max << " km/s, "
            << spectral_grid->wavelength_list.size() << " spectral points\n";
}


VelocityBroadening::~VelocityBroadening()
{
  freeDeviceBuffers();
}


// Rotational kernel (Gray 2005), then a Gaussian of width sigma_kms, each as a
// normalised convolution in velocity space on the actual grid steps.  A kernel
// narrower than every grid step is skipped, since it would leave the spectrum
// unchanged anyway.
void VelocityBroadening::convolveSpectrumCPU(
  const std::vector<double>& spectrum_in,
  std::vector<double>& spectrum_out,
  double sigma_kms,
  double vsini_kms,
  double epsilon)
{
  const bool do_rotation = vsini_kms >= dv_min_kms;
  const bool do_gaussian = 5.0 * sigma_kms >= dv_min_kms;

  if (do_rotation && do_gaussian)
  {
    std::vector<double> temp;
    convolveVelocitySpace(spectrum_in, temp, dv_kms, vsini_kms,
                          RotationalProfile(vsini_kms, epsilon));
    convolveVelocitySpace(temp, spectrum_out, dv_kms, 5.0 * sigma_kms,
                          GaussianProfile(sigma_kms));
  }
  else if (do_rotation)
  {
    convolveVelocitySpace(spectrum_in, spectrum_out, dv_kms, vsini_kms,
                          RotationalProfile(vsini_kms, epsilon));
  }
  else if (do_gaussian)
  {
    convolveVelocitySpace(spectrum_in, spectrum_out, dv_kms, 5.0 * sigma_kms,
                          GaussianProfile(sigma_kms));
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
  if (vsini_kms < dv_min_kms && 5.0 * sigma_total < dv_min_kms)
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
  if (vsini_kms < dv_min_kms && 5.0 * sigma_total < dv_min_kms)
    return;

  const int n_pixels = static_cast<int>(spectral_grid->nbSpectralPoints());

  if (temp_buffer_gpu == nullptr)
    allocateOnDevice(temp_buffer_gpu, static_cast<size_t>(n_pixels));

  if (temp_buffer2_gpu == nullptr)
    allocateOnDevice(temp_buffer2_gpu, static_cast<size_t>(n_pixels));

  if (dv_gpu == nullptr)
  {
    std::vector<float> dv_f(dv_kms.begin(), dv_kms.end());
    moveToDevice(dv_gpu, dv_f);
  }

  // Copy input to temp buffer, then convolve temp -> spectrum_gpu (in-place)
  copyOnDevice(temp_buffer_gpu, spectrum_gpu, static_cast<size_t>(n_pixels));

  applyHighResConvolutionGPU(
    temp_buffer_gpu,
    spectrum_gpu,
    n_pixels,
    dv_gpu,
    dv_min_kms,
    sigma_total,
    vsini_kms,
    epsilon,
    temp_buffer2_gpu);
}


}
