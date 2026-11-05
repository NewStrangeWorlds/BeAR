// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

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

