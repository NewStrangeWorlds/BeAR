/*
* This file is part of the BeAR code (https://github.com/exoclime/BeAR).
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


#ifndef _stellar_activity_h
#define _stellar_activity_h

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

#include "module.h"

#include "../../spectral_grid/spectral_grid.h"
#include "../atmosphere/atmosphere.h"
#include "../stellar_spectrum/stellar_spectrum.h"
#include "../../CUDA_kernels/data_management_kernels.h"

namespace helios {


class StellarContamination : public Module{
  public:
    StellarContamination (
      const std::vector<std::string>& stellar_model_parameters,
      SpectralGrid* spectral_grid_);
    virtual ~StellarContamination() {
      deleteFromDevice(spectrum_phot_gpu);
      deleteFromDevice(spectrum_fac_gpu);
      deleteFromDevice(spectrum_spot_gpu);
    }
    
    virtual void modifySpectrum(
      const std::vector<double>& parameter,
      Atmosphere* atmosphere,
      std::vector<double>& spectrum);
    
    virtual void modifySpectrumGPU(
      const std::vector<double>& parameter,
      Atmosphere* atmosphere,
      double* spectrum_gpu);
  protected:
    SpectralGrid* spectral_grid;
    StellarSpectrumModel* stellar_model;

    double* spectrum_phot_gpu = nullptr;
    double* spectrum_fac_gpu = nullptr;
    double* spectrum_spot_gpu = nullptr;

    size_t nb_stellar_model_param = 0;
};


}

#endif