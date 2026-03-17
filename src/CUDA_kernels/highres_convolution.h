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


#ifndef HIGHRES_CONVOLUTION_H
#define HIGHRES_CONVOLUTION_H


namespace bear {


void applyHighResConvolutionGPU(
  float* spectrum_in_dev,
  float* spectrum_out_dev,
  int    n_pixels,
  double sigma_kms,
  double vsini_kms,
  double delta_v_kms,
  double epsilon,
  float* temp_dev = nullptr);


}


#endif
