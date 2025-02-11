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


#ifndef _spectral_band_h
#define _spectral_band_h

#include <vector>
#include <string>
#include <iostream>

#include "spectral_band_type.h"


namespace bear {

//forward declaration
class SpectralGrid;
class GlobalConfig;



class SpectralBands{
  public:
    SpectralBands(
      GlobalConfig* config_,
      SpectralGrid* spectral_grid_)
      : config(config_),
        spectral_grid(spectral_grid_)
      {}
    ~SpectralBands();
    band_type::id bandType() const {return band_type;}
    void init(
      const std::vector<double>& obs_wavelength_range,
      const std::vector< std::vector<double> >& band_edges,
      const std::vector<double>& band_centres,
      const band_type::id type);
    void init();

    void initDeviceMemory();
    void setInstrumentProfileFWHW(std::vector<double>& profile_fwhm);

    void setBandEdgeIndices(std::vector<double>& wavenumber_grid);
    
    std::vector<double> bandIntegrateSpectrum(
      const std::vector<double>& spectrum, 
      const bool is_flux,
      const bool use_filter_transmission);
    void bandIntegrateSpectrumGPU(
      double* spectrum, 
      double* spectrum_bands, 
      const bool is_flux,
      const bool use_filter_transmission);

    std::vector<double> convolveSpectrum(const std::vector<double>& spectrum);
    void convolveSpectrumGPU(double* spectrum, double* spectrum_processed_dev);

    size_t nbBands() {return nb_bands;}

    std::pair<double, double> obs_wavelength_range = {0.0, 0.0};          //the wavelength range required for the observation
    std::pair<double, double> obs_wavenumber_range = {0.0, 0.0};          //the wavenumber range required for the observation
    std::pair<size_t , size_t> obs_index_range = {0, 0};                  //the spectral indices of the observational range
    
    std::vector<double> center_wavelengths;                               //center wavelengths for each spectral bin
    std::vector< std::vector<double> > edge_wavelengths;                  //wavelengths of the spectral bin edges
    std::vector< std::vector<double> > edge_wavenumbers;                  //wavelengths of the spectral bin edges
    std::vector<std::vector<size_t>> edge_indices;                        //spectral indices corresponding to the edges of the spectral bins

    std::vector<double> instrument_profile_sigma;                         //standard deviation of the instruments profile
    std::vector< std::vector<size_t> > convolution_quadrature_intervals;  //pre-determined limits (indices) for the convolution integration
    
    int* edge_indices_dev;
    int* convolution_start_dev = nullptr;
    int* convolution_end_dev = nullptr;
    double* instrument_profile_sigma_dev = nullptr;
    int* band_start_dev = nullptr;
    int* band_end_dev = nullptr;

  private:
    GlobalConfig* config = nullptr;
    SpectralGrid* spectral_grid = nullptr;
    band_type::id band_type;
    size_t nb_bands = 0;                                                   //number of sub-bands/bins 

    double bandIntegrateSpectrumFlux(
      const std::vector<double>& spectrum,
      const size_t& band);
    double bandIntegrateSpectrum(
      const std::vector<double>& spectrum, 
      const size_t& band,
      const bool use_filter_transmission);
    double convolveSpectrum(
      const std::vector<double>& spectrum,
      const unsigned int index);
    void setConvolutionQuadratureIntervals();
    void setConvolutionQuadratureIntervals(
      const size_t index,
      const double cutoff_distance);
};




}


#endif
