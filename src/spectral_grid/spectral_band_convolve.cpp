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
#include <algorithm>
#include <cmath>


#include "spectral_band.h"
#include "spectral_grid.h"
#include "../additional/quadrature.h"
#include "../additional/aux_functions.h"


namespace helios{


//Since we want to do the convolution only with a certain distance from the centre, we here
//pre-determine the edges of the quadruature intervalls for each wavenumber
//we use a distance of 5 sigma here
void SpectralBands::setConvolutionQuadratureIntervals()
{
  const size_t nb_high_res_points = obs_index_range.second - obs_index_range.first + 1;

  convolution_quadrature_intervals.assign(nb_high_res_points, std::vector<size_t>(2,0));

  for (size_t i=0; i<nb_high_res_points; ++i)
  {
    //initialise with largest possible interval
    convolution_quadrature_intervals[i][0] = obs_index_range.first;
    convolution_quadrature_intervals[i][1] = obs_index_range.second;

    const double cutoff_distance = 5.0 * instrument_profile_sigma[i];

    setConvolutionQuadratureIntervals(i+obs_index_range.first, cutoff_distance);
  }
}



//find the edges of the quadruature intervals for a wavenumber at a certain distance
//the indices refer to the position in the global spectral grid
void SpectralBands::setConvolutionQuadratureIntervals(
  const size_t index,
  const double cutoff_distance)
{ 
  const size_t conv_intervall_idx = index - obs_index_range.first;

  if (cutoff_distance == 0)
  {
    convolution_quadrature_intervals[conv_intervall_idx][0] = index;
    convolution_quadrature_intervals[conv_intervall_idx][1] = index;
  }
  else
  {
    const double left_wavelength_edge = spectral_grid->wavelength_list[index] + cutoff_distance;
    convolution_quadrature_intervals[conv_intervall_idx][0] = spectral_grid->findClosestIndex(
      left_wavelength_edge,
      spectral_grid->wavelength_list,
      spectral_grid->wavelength_list.begin());

    const double right_wavelength_edge = spectral_grid->wavelength_list[index] - cutoff_distance;
    convolution_quadrature_intervals[conv_intervall_idx][1] = spectral_grid->findClosestIndex(
      right_wavelength_edge,
      spectral_grid->wavelength_list,
      spectral_grid->wavelength_list.begin() + convolution_quadrature_intervals[conv_intervall_idx][0]);
  }

}



//convolve the spectrum with an instrument profile
std::vector<double> SpectralBands::convolveSpectrum(const std::vector<double>& spectrum)
{
  if (instrument_profile_sigma.size() == 0) return spectrum;

  std::vector<double> convolved_spectrum(spectrum.size(), 0.0);

  #pragma omp parallel for
  for (size_t i=obs_index_range.first; i<obs_index_range.second+1; ++i)
    convolved_spectrum[i] = convolveSpectrum(spectrum, i);


  return convolved_spectrum;
}



//calculate the convolved spectrum at a given spectral index
double SpectralBands::convolveSpectrum(
  const std::vector<double>& spectrum,
  const unsigned int index)
{
  const size_t conv_idx = index - obs_index_range.first;

  if (convolution_quadrature_intervals[conv_idx][0] == convolution_quadrature_intervals[conv_idx][1])
    return spectrum[index];

  //copy the relevant part from the wavelength vector
  const std::vector<double> x(
    spectral_grid->wavelength_list.begin()+convolution_quadrature_intervals[conv_idx][0],
    spectral_grid->wavelength_list.begin()+convolution_quadrature_intervals[conv_idx][1]);

  std::vector<double> y(
    spectrum.begin()+convolution_quadrature_intervals[conv_idx][0],
    spectrum.begin()+convolution_quadrature_intervals[conv_idx][1]);


  for (size_t i=0; i<y.size(); ++i)
  {
    const double distance = std::abs(spectral_grid->wavelength_list[index] - x[i]);

    y[i] *= aux::normalDistribution(instrument_profile_sigma[index], distance);
  }


  return std::abs(aux::quadratureTrapezoidal(x, y));
}



}
