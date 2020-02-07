/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2020 Daniel Kitzmann
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



namespace helios{



void SpectralBands::init(GlobalConfig* global_config, SpectralGrid* grid, 
                        const std::vector< std::vector<double> >& band_edges, const std::vector<double>& band_centres, 
                        const BandType type, const std::string filter_file_name)
{
  config = global_config;
  spectral_grid = grid;

  band_type = type;
  
  spectral_grid->sampleWavelengths(band_edges, global_config->spectral_resolution, global_spectral_indices, edge_indices);

  nb_bands = edge_indices.size() - 1;

 
  band_centers_wavelength = band_centres;
  band_edges_wavelength = band_edges;
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




void SpectralBands::setInstrumentProfileFWHW(std::vector<double>& profile_fwhm)
{
  if (profile_fwhm.size() == 0)
    return;


  if (profile_fwhm.size() != band_centers_wavelength.size())
  {
    std::cout << "Profile FWHM vector has incorrect size \n";
    return;
  }


  instrument_profile_sigma.assign(wavenumbers.size(), 0.0);
  

  for (size_t i=0; i<wavenumbers.size(); ++i)
  {
    for (size_t j=0; j<nb_bands; ++j)
    {
      if (wavenumbers[i] >= band_wavenumbers[j].front() && wavenumbers[i] <= band_wavenumbers[j].back())
        instrument_profile_sigma[i] = profile_fwhm[j] / 2.355;  //convert FWHM to standard deviation
    }
  }


  setConvolutionQuadratureIntervals();
}



}
