// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _stellar_activity_h
#define _stellar_activity_h

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <memory>

#include "../module.h"

#include "../../../spectral_grid/spectral_grid.h"
#include "../../atmosphere/atmosphere.h"
#include "../../stellar_spectrum/stellar_spectrum.h"
#include "../../../CUDA_kernels/data_management_kernels.h"

namespace bear {


class StellarContamination : public Module{
  public:
    StellarContamination (
      const std::vector<std::string>& stellar_model_parameters,
      SpectralGrid* spectral_grid_);
    virtual ~StellarContamination() {
      if (spectrum_phot_gpu != nullptr) deleteFromDevice(spectrum_phot_gpu);
      if (spectrum_fac_gpu != nullptr) deleteFromDevice(spectrum_fac_gpu);
      if (spectrum_spot_gpu != nullptr) deleteFromDevice(spectrum_spot_gpu);
    }
    
    virtual void modifySpectrum(
      const std::vector<double>& parameter,
      Atmosphere* atmosphere,
      std::vector<double>& spectrum);
    
    virtual void modifySpectrumGPU(
      const std::vector<double>& parameter,
      Atmosphere* atmosphere,
      float* spectrum_gpu);
  protected:
    SpectralGrid* spectral_grid;
    std::unique_ptr<StellarSpectrumModel> stellar_model;

    float* spectrum_phot_gpu = nullptr;
    float* spectrum_fac_gpu = nullptr;
    float* spectrum_spot_gpu = nullptr;

    size_t nb_stellar_model_param = 0;
};


}

#endif