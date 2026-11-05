// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-FileCopyrightText: 2026 Esther Wang
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _velocity_broadening_h
#define _velocity_broadening_h

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <memory>

#include "../module.h"

#include "../../../spectral_grid/spectral_grid.h"
#include "../../../CUDA_kernels/data_management_kernels.h"

namespace bear {


class VelocityBroadening : public Module{
  public:
    VelocityBroadening (
      const std::vector<std::string>& velocity_broadening_parameters,
      SpectralGrid* spectral_grid_);
    virtual ~VelocityBroadening();

    virtual void modifySpectrum(
      const std::vector<double>& parameter,
      Atmosphere* atmosphere,
      std::vector<double>& spectrum);

    virtual void modifySpectrumGPU(
      const std::vector<double>& parameter,
      Atmosphere* atmosphere,
      float* spectrum_gpu);

    void setSpectralGrid(SpectralGrid* grid);
  protected:
    SpectralGrid* spectral_grid;

    // Velocity step between neighbouring grid points (km/s), and its minimum.
    // The high-res grid is not uniform in ln(lambda), so the kernels use these
    // actual steps instead of a single value for the whole grid.
    std::vector<double> dv_kms;
    double dv_min_kms = 0;
    float* dv_gpu = nullptr;

    float* temp_buffer_gpu = nullptr;
    float* temp_buffer2_gpu = nullptr;

    void setVelocitySteps();
    void freeDeviceBuffers();

    void convolveSpectrumCPU(
      const std::vector<double>& spectrum_in,
      std::vector<double>& spectrum_out,
      double sigma_kms,
      double vsini_kms,
      double epsilon);
};


}

#endif