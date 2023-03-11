/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2022 Daniel Kitzmann
*
* Helios-r2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Helios-r2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* Helios-r2 directory under <LICENSE>. If not, see
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


namespace helios{


void SpectralBands::init(
  const std::vector<double>& extended_edges, 
  const std::vector< std::vector<double> >& obs_band_edges, 
  const std::vector<double>& obs_band_centres, 
  const BandType type)
{
  band_type = type;

  if (extended_edges.size() > 0)
  {
    std::vector< std::vector<double> > band_edges;
    band_edges.reserve(obs_band_edges.size() + 2);

    band_edges.push_back(std::vector<double> {extended_edges[0], obs_band_edges[0][0]});
    band_edges.insert(band_edges.end(), obs_band_edges.begin(), obs_band_edges.end());
    band_edges.push_back(std::vector<double> {obs_band_edges.back()[1], extended_edges[1]});
    
    std::vector<size_t> all_edge_indices;

    spectral_grid->sampleWavelengths(band_edges, config->spectral_resolution, global_spectral_indices, all_edge_indices);

    edge_indices = {all_edge_indices.begin()+1, all_edge_indices.end() - 1};
  }
  else
    spectral_grid->sampleWavelengths(obs_band_edges, config->spectral_resolution, global_spectral_indices, edge_indices);

  nb_bands = edge_indices.size() - 1;

 
  band_centers_wavelength = obs_band_centres;
  band_edges_wavelength = obs_band_edges;
}


//retrieves the local indices and the wavenumbers for this band
//only to be called once the full spectral grid has been assembled!
void SpectralBands::setLocalIndices()
{
  spectral_indices = spectral_grid->globalToLocalIndex(global_spectral_indices);
  wavenumbers = spectral_grid->wavenumberList(spectral_indices);

  spectral_grid->convertWavenumbersToWavelengths(wavenumbers, wavelengths);


  band_spectral_indices.resize(nb_bands);
  band_wavenumbers.resize(nb_bands);


  //setting up the structures for the sub-bands
  for (size_t i=0; i<nb_bands; ++i)
  {
    for (size_t j=edge_indices[i]; j<edge_indices[i+1]+1; ++j )
      band_spectral_indices[i].push_back(spectral_indices[j]);

    band_wavenumbers[i] = spectral_grid->wavenumberList(band_spectral_indices[i]);
  }

}




void SpectralBands::setInstrumentProfileFWHW(std::vector<double>& obs_profile_fwhm)
{
  if (obs_profile_fwhm.size() == 0)
    return;


  if (obs_profile_fwhm.size() != band_centers_wavelength.size())
  {
    std::cout << "Profile FWHM vector has incorrect size \n";
    return;
  }
  
  std::vector<double> profile_fwhm = spectral_grid->interpolateToWavelengthGrid(
    band_centers_wavelength,
    obs_profile_fwhm,
    wavelengths,
    false);


  instrument_profile_sigma.assign(wavenumbers.size(), 0.0);


  for (size_t i=0; i<wavelengths.size(); ++i)
    instrument_profile_sigma[i] = profile_fwhm[i]/ 2.355;  //convert FWHM to standard deviation


  setConvolutionQuadratureIntervals();
}



void SpectralBands::initDeviceMemory()
{
  std::vector<int> band_indices(spectral_indices.begin(), spectral_indices.end());
  moveToDevice(spectral_indices_dev, band_indices);

  moveToDevice(wavelengths_dev, wavelengths);

  std::vector<int> band_start(nb_bands, 0);
  std::vector<int> band_end(nb_bands, 0);

  for (size_t i=0; i<nb_bands; ++i)
  {
    for (size_t j=0; j<spectral_indices.size(); ++j)
      if (band_spectral_indices[i].front() == spectral_indices[j])
      {
        band_start[i] = j;
        break;
      }

    for (size_t j=0; j<spectral_indices.size(); ++j)
      if (band_spectral_indices[i].back() == spectral_indices[j])
      {
        band_end[i] = j;
        break;
      }
  }

  moveToDevice(band_start_dev, band_start);
  moveToDevice(band_end_dev, band_end);


  if (instrument_profile_sigma.size() > 0)
  {
    std::vector<int> convolution_start(wavelengths.size(), 0);
    std::vector<int> convolution_end(wavelengths.size(), 0);

    for (size_t i=0; i<wavelengths.size(); ++i)
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
  deleteFromDevice(spectral_indices_dev);
  deleteFromDevice(wavelengths_dev);

  deleteFromDevice(band_start_dev);
  deleteFromDevice(band_end_dev);


  deleteFromDevice(convolution_start_dev);
  deleteFromDevice(convolution_end_dev);
  deleteFromDevice(instrument_profile_sigma_dev);
}



}
