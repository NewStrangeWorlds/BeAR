// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _star_spectrum_file_h
#define _star_spectrum_file_h

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

#include "stellar_spectrum.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../CUDA_kernels/data_management_kernels.h"


namespace bear {


class StarSpectrumFile : public StellarSpectrumModel{
  public:
    StarSpectrumFile (
      const std::string file_path,
      SpectralGrid* spectral_grid_);
    StarSpectrumFile (
      const std::vector<double>& stellar_spectrum_wavelengths,
      const std::vector<double>& stellar_spectrum_flux,
      SpectralGrid* spectral_grid_);
    virtual ~StarSpectrumFile() {
      if (spectrum_dev != nullptr) deleteFromDevice(spectrum_dev);
    }
    
    virtual std::vector<double> calcFlux(
      const std::vector<double>& parameter);
    
    virtual void calcFluxGPU(
      const std::vector<double>& parameter,
      float* spectrum_gpu);
  protected:
    SpectralGrid* spectral_grid;
    std::vector<double> spectrum;
    float* spectrum_dev = nullptr;

    void readSpectrum(
      const std::string file_path,
      std::vector<double>& spectrum_file,
      std::vector<double>& wavelength_file);
};


}

#endif