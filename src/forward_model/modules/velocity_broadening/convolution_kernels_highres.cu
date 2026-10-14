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
 * convolution_kernels_highres.cu
 *
 * GPU spectral broadening for high-resolution spectroscopy.
 * Two-pass pipeline: rotational broadening (Gray 2005) then instrumental
 * Gaussian broadening, both as normalised convolutions in velocity space that
 * use the actual velocity step between neighbouring grid points (see
 * velocity_convolution_kernels.cuh).
 */


#include <cmath>
#include <stdio.h>

#include "../../../CUDA_kernels/error_check.h"
#include "../../../CUDA_kernels/data_management_kernels.h"
#include "velocity_convolution_kernels.cuh"
#include "../../../CUDA_kernels/highres_convolution.h"


namespace bear {


__host__
void applyHighResConvolutionGPU(
  float*       spectrum_in_dev,
  float*       spectrum_out_dev,
  int          n_pixels,
  const float* dv_dev,
  double       dv_min_kms,
  double       sigma_kms,
  double       vsini_kms,
  double       epsilon,
  float*       temp_dev)
{
  cudaGetLastError();

  // A kernel whose support is narrower than every step of the grid leaves the
  // spectrum unchanged, so these tests are exact rather than approximations.
  const bool do_rotation = vsini_kms >= dv_min_kms;
  const bool do_gaussian = 5.0 * sigma_kms >= dv_min_kms;

  const RotationalProfileGPU rotation{
    static_cast<float>(1.0 / vsini_kms),
    static_cast<float>(2.0 * (1.0 - epsilon)),
    static_cast<float>(0.5 * M_PI * epsilon)};

  const GaussianProfileGPU gaussian{
    static_cast<float>(1.0 / (2.0 * sigma_kms * sigma_kms))};

  const float rot_hw   = static_cast<float>(vsini_kms);
  const float gauss_hw = static_cast<float>(5.0 * sigma_kms);

  bool allocated_temp = false;

  if (do_rotation && do_gaussian)
  {
    if (temp_dev == nullptr)
    {
      allocateOnDevice(temp_dev, (size_t)n_pixels);
      allocated_temp = true;
    }

    launchVelocitySpaceConvolution(
      spectrum_in_dev, temp_dev, n_pixels, dv_dev, rot_hw, rotation);
    launchVelocitySpaceConvolution(
      temp_dev, spectrum_out_dev, n_pixels, dv_dev, gauss_hw, gaussian);
  }
  else if (do_rotation)
  {
    launchVelocitySpaceConvolution(
      spectrum_in_dev, spectrum_out_dev, n_pixels, dv_dev, rot_hw, rotation);
  }
  else if (do_gaussian)
  {
    launchVelocitySpaceConvolution(
      spectrum_in_dev, spectrum_out_dev, n_pixels, dv_dev, gauss_hw, gaussian);
  }
  else
  {
    // No broadening — copy input to output
    copyOnDevice(spectrum_out_dev, spectrum_in_dev, (size_t)n_pixels);
    CUDA_CHECK_AFTER_KERNEL();
  }

  if (allocated_temp)
    deleteFromDevice(temp_dev);
}


} // namespace bear
