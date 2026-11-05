// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

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

    float* spectrum_gpu = nullptr;
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
      float* spectrum_gpu);
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