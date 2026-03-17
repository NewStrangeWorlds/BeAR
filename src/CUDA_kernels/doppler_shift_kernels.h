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


#ifndef DOPPLER_SHIFT_KERNELS_H
#define DOPPLER_SHIFT_KERNELS_H


namespace bear {


// Doppler-shift a spectrum on the GPU by linear interpolation.
// On a constant-resolution (log-lambda) grid, a radial velocity shift
// corresponds to a constant pixel shift: shift_pixels = v_rad / delta_v_kms.
void dopplerShiftGPU(
  const float* spectrum_in_dev,
  float* spectrum_out_dev,
  int n_pixels,
  float shift_pixels);

// Interpolate a model spectrum onto a target wavelength grid (GPU).
// Both wavelength arrays must be sorted in ascending order.
void interpolateSpectrumGPU(
  const float* model_spectrum_dev,
  const float* model_wavelengths_dev,
  int nb_model_points,
  const float* target_wavelengths_dev,
  float* target_spectrum_dev,
  int nb_target_points);


}


#endif
