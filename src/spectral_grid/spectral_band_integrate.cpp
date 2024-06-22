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
#include <algorithm>


#include "spectral_grid.h"
#include "../additional/quadrature.h"


namespace bear{


//Integration of the high-res spectrum in each observational band
std::vector<double> SpectralBands::bandIntegrateSpectrum(
  const std::vector<double>& spectrum, 
  const bool is_flux,
  const bool use_filter_transmission)
{
  std::vector<double> band_values(nb_bands, 0.0);

  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t i=0; i<nb_bands; ++i)
  {
    if (is_flux)
      band_values[i] = bandIntegrateSpectrumFlux(
        spectrum, 
        i);
    else
      band_values[i] = bandIntegrateSpectrum(
        spectrum, 
        i,
        use_filter_transmission);
  }

  return band_values;
}


//Integration of the high-res spectrum for a specific band
//Note that the units of the spectrum are W m-2 cm, the integrated band values are in W m-2 mu-1
double SpectralBands::bandIntegrateSpectrumFlux(
  const std::vector<double>& spectrum, const size_t& band)
{
  double wavelength_edge_left =  1.0/spectral_grid->wavenumber_list[edge_indices[band][0]] * 1e4;
  double wavelength_edge_right =  1.0/spectral_grid->wavenumber_list[edge_indices[band][1]] * 1e4;

  double band_mean = 
    aux::quadratureTrapezoidal(
      spectral_grid->wavenumber_list, 
      spectrum,
      edge_indices[band][0],
      edge_indices[band][1]) / (wavelength_edge_left - wavelength_edge_right);

  return band_mean;
}


//Integration of the high-res spectrum for a specific band
//Note that the units of the spectrum are W m-2 cm, the integrated band values are in W m-2 mu-1
double SpectralBands::bandIntegrateSpectrum(
  const std::vector<double>& spectrum, 
  const size_t& band,
  const bool use_filter_transmission)
{
  double wavelength_edge_right =  spectral_grid->wavelength_list[edge_indices[band][1]];
  double wavelength_edge_left = spectral_grid->wavelength_list[edge_indices[band][0]];


  double band_mean = 
    aux::quadratureTrapezoidal(
      spectral_grid->wavelength_list, 
      spectrum,
      edge_indices[band][0],
      edge_indices[band][1]);

  band_mean = std::abs(band_mean);
      
  if (!use_filter_transmission) 
    band_mean /= (wavelength_edge_left - wavelength_edge_right);

  return band_mean;
}



}

