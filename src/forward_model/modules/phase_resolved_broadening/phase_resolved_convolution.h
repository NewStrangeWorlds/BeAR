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
