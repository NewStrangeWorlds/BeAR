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


#ifndef _stellar_spectrum_grid_h
#define _stellar_spectrum_grid_h

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

#include "stellar_spectrum.h"
#include "../../spectral_grid/spectral_grid.h"


namespace bear {


class SpectrumFile {
  public:
    SpectrumFile(
      const std::string file_path_)
        : file_path(file_path_) 
        {}
    void loadFile();
    void unloadData();
    
    const std::string file_path = "";
    bool is_loaded = false;

    std::vector<double> spectrum;
};


class SampledStellarSpectrum{
  public:
    SampledStellarSpectrum(
      const double temperature_, 
      const double log_g_,
      const double metallicity_,
      const std::string& file_path)
        : temperature(temperature_)
        , log_temperature(std::log10(temperature_))
        , log_g(log_g_)
        , metallicity(metallicity_)
        , data_file(file_path) 
        {}
    ~SampledStellarSpectrum();
    void sampleSpectrum(
      SpectralGrid* spectral_grid,
      const std::vector<double>& grid_wavelengths,
      const bool use_gpu);
    void deleteSampledData();
    
    const double temperature = 0.0;
    const double log_temperature = 0.0;
    const double log_g = 0.0;
    const double metallicity = 0.0;
    
    bool is_sampled = false;

    double* spectrum_gpu = nullptr;
    std::vector<double> spectrum;
  private:
    SpectrumFile data_file;
};



class StellarSpectrumGrid : public StellarSpectrumModel{
  public:
    StellarSpectrumGrid (
      const std::string file_path,
      SpectralGrid* spectral_grid_);
    virtual ~StellarSpectrumGrid() {}
    
    virtual std::vector<double> calcFlux(
      const std::vector<double>& parameter);
    
    virtual void calcFluxGPU(
      const std::vector<double>& parameter,
      double* spectrum_gpu);
  protected:
    SpectralGrid* spectral_grid;
    std::vector<double> grid_wavelengths;
    std::string wavelength_grid_file_name = "";
    bool is_wavelength_grid_reversed = false;

    std::vector<double> effective_temperature;
    std::vector<double> log_g;
    std::vector<double> metallicity;

    std::vector<SampledStellarSpectrum> grid_spectra;
    std::vector<std::vector<std::vector<SampledStellarSpectrum*>>> grid; 

    std::vector<std::string> loadParameterFile(const std::string& file_path);
    void loadWavelengthGrid(const std::string& file_path);
    void parseFileList(
      const std::vector<std::string>& file_list,
      const std::string& folder);
    
    std::pair<size_t, size_t> findClosestGridPoints(
      const std::vector<double>& grid_points, 
      const double point);
};


}

#endif