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


#include "spectral_band.h"

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "spectral_band_type.h"
#include "spectral_grid.h"
#include "../config/global_config.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../additional/exceptions.h"


namespace bear{


void SpectralBands::init(
  const std::vector<double>& obs_wavelength_range_, 
  const std::vector< std::vector<double> >& obs_band_edges, 
  const std::vector<double>& obs_band_centres, 
  const BandType type)
{
  band_type = type;
  center_wavelengths = obs_band_centres;

  edge_wavelengths = obs_band_edges;
  edge_wavenumbers = edge_wavelengths;

  for (auto & w : edge_wavenumbers)
  {
    w[0] = spectral_grid->wavelengthToWavenumber(w[0]);
    w[1] = spectral_grid->wavelengthToWavenumber(w[1]);
  }

  obs_wavelength_range.first = obs_wavelength_range_[0];
  obs_wavelength_range.second = obs_wavelength_range_[1];

  obs_wavenumber_range.first = spectral_grid->wavelengthToWavenumber(obs_wavelength_range.first);
  obs_wavenumber_range.second = spectral_grid->wavelengthToWavenumber(obs_wavelength_range.second);

  nb_bands = center_wavelengths.size();
}


void SpectralBands::init()
{
  obs_index_range.first = spectral_grid->findClosestIndex(
    obs_wavenumber_range.first,
    spectral_grid->wavenumber_list,
    spectral_grid->wavenumber_list.begin());
  obs_index_range.second = spectral_grid->findClosestIndex(
    obs_wavenumber_range.second,
    spectral_grid->wavenumber_list,
    spectral_grid->wavenumber_list.begin());
}



void SpectralBands::setBandEdgeIndices(std::vector<double>& wavenumber_grid)
{
  edge_indices.resize(0);
  edge_indices.reserve(nb_bands);

  auto it_start = wavenumber_grid.begin();
  double last_wavenumber = 0;

  for (auto & b : edge_wavenumbers)
  {
    //in case of overlapping bands, we need to start the search at the beginning
    if (last_wavenumber > b[0]) 
      it_start = wavenumber_grid.begin();
    
    const size_t idx_1 = spectral_grid->findClosestIndex(b[0], wavenumber_grid, it_start);
    const size_t idx_2 = spectral_grid->findClosestIndex(b[1], wavenumber_grid, it_start);
      
    if (idx_1 > wavenumber_grid.size() -1 || idx_2 > wavenumber_grid.size() -1)
    {
      std::string error_message = "Edge wavenumber index not found!\t" 
        + std::to_string(b[0]) + "  " + std::to_string(b[1]) + "\n";

      throw InvalidInput(std::string ("SpectralBand::setBandEdgeIndices"), error_message);
    }
    
    if (idx_1 == idx_2)
    {
      std::string error_message = "The spectral resolution of the opacity grid is too small!\nUnable to locate enough points to resolve the bin from " 
        + std::to_string(1.0/b[0]*1e4) + " microns to " + std::to_string(1.0/b[1]*1e4) + " microns\n";

      throw InvalidInput(std::string ("SpectralBand::setBandEdgeIndices"), error_message);

    }
    edge_indices.push_back(std::vector<size_t>{idx_1, idx_2});

    last_wavenumber = b[1];
    it_start = wavenumber_grid.begin() + idx_2;
  }

  init();
  //std::cout << "\n";
}


//interpolates the instrument profile onto the global wavelength grid
//the instrument profile has the same size as the global grid
//but will be 0 outside of the observational range
void SpectralBands::setInstrumentProfileFWHW(std::vector<double>& obs_profile_fwhm)
{
  if (obs_profile_fwhm.size() == 0)
    return;

  if (obs_profile_fwhm.size() != center_wavelengths.size())
  {
    std::cout << "Profile FWHM vector has incorrect size \n";
    return;
  }

  std::vector<double> high_res_wavelengths(
    spectral_grid->wavelength_list.begin()+obs_index_range.first,
    spectral_grid->wavelength_list.begin()+obs_index_range.second+1);

  std::vector<double> profile_fwhm = spectral_grid->interpolateToWavelengthGrid(
    center_wavelengths,
    obs_profile_fwhm,
    high_res_wavelengths,
    false);

  instrument_profile_sigma.assign(spectral_grid->nbSpectralPoints(), 0.0);

  for (size_t i=0; i<high_res_wavelengths.size(); ++i)
    instrument_profile_sigma[i+obs_index_range.first] = profile_fwhm[i]/ 2.355;  //convert FWHM to standard deviation

  setConvolutionQuadratureIntervals();
}



void SpectralBands::initDeviceMemory()
{
  if (config->use_gpu == false) 
    return;

  std::vector<int> band_start(nb_bands, 0);
  std::vector<int> band_end(nb_bands, 0);

  for (size_t i=0; i<nb_bands; ++i)
  {
    band_start[i] = edge_indices[i][0];
    band_end[i] = edge_indices[i][1];
  }

  moveToDevice(band_start_dev, band_start);
  moveToDevice(band_end_dev, band_end);

  const size_t nb_high_res_points = obs_index_range.second - obs_index_range.first + 1;

  if (instrument_profile_sigma.size() > 0)
  {
    std::vector<int> convolution_start(nb_high_res_points, 0);
    std::vector<int> convolution_end(nb_high_res_points, 0);

    for (size_t i=0; i<nb_high_res_points; ++i)
    {
      convolution_start[i] = convolution_quadrature_intervals[i][0];
      convolution_end[i] = convolution_quadrature_intervals[i][1];
    }

    moveToDevice(convolution_start_dev, convolution_start);
    moveToDevice(convolution_end_dev, convolution_end);

    moveToDevice(instrument_profile_sigma_dev, instrument_profile_sigma);
  }

}


SpectralBands::~SpectralBands()
{
  if (config->use_gpu)
  {
    deleteFromDevice(band_start_dev);
    deleteFromDevice(band_end_dev);

    deleteFromDevice(convolution_start_dev);
    deleteFromDevice(convolution_end_dev);
    deleteFromDevice(instrument_profile_sigma_dev);
  }

}



}
