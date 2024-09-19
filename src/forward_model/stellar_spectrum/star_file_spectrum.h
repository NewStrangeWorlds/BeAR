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
      deleteFromDevice(spectrum_dev);
    }
    
    virtual std::vector<double> calcFlux(
      const std::vector<double>& parameter);
    
    virtual void calcFluxGPU(
      const std::vector<double>& parameter,
      double* spectrum_gpu);
  protected:
    SpectralGrid* spectral_grid;
    std::vector<double> spectrum;
    double* spectrum_dev = nullptr;

    void readSpectrum(
      const std::string file_path,
      std::vector<double>& spectrum_file,
      std::vector<double>& wavelength_file);
};


}

#endif