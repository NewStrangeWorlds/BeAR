/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2019 Daniel Kitzmann
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
#include <algorithm>


#include "spectral_grid.h"
#include "../additional/quadrature.h"


namespace helios{


//Integration of the high-res spectrum in each observational band
void SpectralBands::bandIntegrateSpectrum(const std::vector<double>& spectrum, std::vector<double>& band_values)
{

  band_values.assign(nb_bands, 0.0);

  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t i=0; i<nb_bands; ++i)
    band_values[i] = bandIntegrateSpectrum(spectrum, i);
}


//Integration of the high-res spectrum for a specific band
//Note that the units of the spectrum are W m-2 cm, the integrated band values are in W m-2 mu-1
double SpectralBands::bandIntegrateSpectrum(const std::vector<double>& spectrum, const size_t& band)
{

  std::vector<double> spectrum_subset(band_spectral_indices[band].size(), 0.0);

  for (size_t i=0; i<band_spectral_indices[band].size(); ++i)
    spectrum_subset[i] = spectrum[ band_spectral_indices[band][i] ];


  double wavelength_edge_left =  1.0/band_wavenumbers[band].front() * 1e4;
  double wavelength_edge_right =  1.0/band_wavenumbers[band].back() * 1e4;


  double band_mean = aux::quadratureTrapezoidal(band_wavenumbers[band], spectrum_subset) / (wavelength_edge_left - wavelength_edge_right);


  return band_mean;
}



}

