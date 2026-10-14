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
 * convolution_kernels_phase_resolved.cu
 *
 * GPU convolution with the tabulated broadening kernel of the phase-resolved
 * broadening module (Brogi et al. 2016), as a normalised convolution in
 * velocity space on the actual grid steps (velocity_convolution_kernels.cuh).
 */


#include <cmath>
#include <cstdio>

#include "../../../CUDA_kernels/error_check.h"
#include "../velocity_broadening/velocity_convolution_kernels.cuh"
#include "phase_resolved_convolution.h"


namespace bear {


__host__
void applyTabulatedConvolutionGPU(
  const float* spectrum_in_dev,
  float*       spectrum_out_dev,
  int          n_pixels,
  const float* dv_dev,
  const float* table_dev,
  int          table_size,
  double       table_v_min,
  double       table_step,
  double       half_width)
{
  cudaGetLastError();

  const TabulatedProfileGPU profile{
    table_dev,
    table_size,
    static_cast<float>(table_v_min),
    static_cast<float>(1.0 / table_step)};

  launchVelocitySpaceConvolution(
    spectrum_in_dev, spectrum_out_dev, n_pixels, dv_dev,
    static_cast<float>(half_width), profile);
}


}
