// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-FileCopyrightText: 2026 Esther Wang
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

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
