// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-FileCopyrightText: 2026 Esther Wang
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef HIGHRES_CONVOLUTION_H
#define HIGHRES_CONVOLUTION_H


namespace bear {


// Rotational (Gray 2005) then Gaussian broadening on the high-resolution grid,
// as normalised convolutions in velocity space.  dv_dev holds the n_pixels - 1
// velocity steps between neighbouring grid points (km/s), dv_min_kms the
// smallest of them.
void applyHighResConvolutionGPU(
  float*       spectrum_in_dev,
  float*       spectrum_out_dev,
  int          n_pixels,
  const float* dv_dev,
  double       dv_min_kms,
  double       sigma_kms,
  double       vsini_kms,
  double       epsilon,
  float*       temp_dev = nullptr);


}


#endif
