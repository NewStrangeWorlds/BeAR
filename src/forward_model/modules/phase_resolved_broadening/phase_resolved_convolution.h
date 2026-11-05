// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-FileCopyrightText: 2026 Esther Wang
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef PHASE_RESOLVED_CONVOLUTION_H
#define PHASE_RESOLVED_CONVOLUTION_H


namespace bear {


// Normalised convolution in velocity space with a kernel tabulated at the
// velocities table_v_min + k * table_step (k = 0..table_size-1), evaluated at the
// actual velocity offsets of the grid points (dv_dev: the n_pixels - 1 steps
// between neighbouring points, km/s).
void applyTabulatedConvolutionGPU(
  const float* spectrum_in_dev,
  float*       spectrum_out_dev,
  int          n_pixels,
  const float* dv_dev,
  const float* table_dev,
  int          table_size,
  double       table_v_min,
  double       table_step,
  double       half_width);


}


#endif
