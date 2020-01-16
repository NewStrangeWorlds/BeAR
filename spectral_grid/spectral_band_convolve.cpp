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
#include <cmath>


#include "spectral_band.h"
#include "../additional/quadrature.h"
#include "../additional/aux_functions.h"


namespace helios{


//Since we want to do the convolution only with a certain distance from the centre, we here
//pre-determine the edges of the quadruature intervalls for each wavenumber
//we use a distance of 5 sigma here
void SpectralBands::setConvolutionQuadratureIntervals()
{
  convolution_quadrature_intervals.assign(wavelengths.size(), std::vector<size_t>(2,0));

  
  for (size_t i=0; i<wavelengths.size(); ++i)
  {
    //initialise with largest possible interval
    convolution_quadrature_intervals[i][0] = 0;
    convolution_quadrature_intervals[i][1] = wavelengths.size()-1;

    const double cutoff_distance = 5.0 * instrument_profile_sigma[i];

    setConvolutionQuadratureIntervals(i, cutoff_distance);    
  }
  
}



//find the edges of the quadruature intervalls for a wavenumber at a certain distance
void SpectralBands::setConvolutionQuadratureIntervals(const unsigned int index, const double cutoff_distance)
{

  for (size_t j=0; j< wavelengths.size()-1; ++j)
  {
    if ( (wavelengths[j] - wavelengths[index]) > cutoff_distance && (wavelengths[j+1] - wavelengths[index]) < cutoff_distance )
      convolution_quadrature_intervals[index][0] = j;
  }


  for (size_t j=convolution_quadrature_intervals[index][0]; j< wavelengths.size()-1; ++j)
  {
    if ( (wavelengths[index] - wavelengths[j]) < cutoff_distance && (wavelengths[index] - wavelengths[j+1]) > cutoff_distance )
      convolution_quadrature_intervals[index][1] = j;
  }


}



//convolve the spectrum with an instrument profile
void SpectralBands::convolveSpectrum(const std::vector<double>& spectrum, std::vector<double>& convolved_spectrum)
{
  
  if (convolution_quadrature_intervals.size() == 0)
    setConvolutionQuadratureIntervals();

  const size_t nb_global_wavenumbers = spectrum.size();
  const size_t nb_wavenumbers = wavenumbers.size();
  
  //extract the relevant wavelength range of the band from the full spectrum 
  std::vector<double> band_spectrum(nb_wavenumbers, 0.0);

  for (size_t i=0; i<nb_wavenumbers; ++i) 
    band_spectrum[i] = spectrum[ spectral_indices[i] ];
  
  
  convolved_spectrum.assign(nb_global_wavenumbers, 0.0);

  #pragma omp parallel for
  for (size_t i=0; i<nb_wavenumbers; ++i)
    convolved_spectrum[spectral_indices[i]] = convolveSpectrum(band_spectrum, i);

 
}



//calculate the convolved spectrum at a given spectral index
double SpectralBands::convolveSpectrum(const std::vector<double>& spectrum, const unsigned int index)
{
  const size_t lower = convolution_quadrature_intervals[index][0];
  const size_t upper = convolution_quadrature_intervals[index][1];
  const size_t length = (upper - lower) + 1;


  //copy the relevant part from the wavelength vector
  const std::vector<double> x(wavelengths.begin()+lower, wavelengths.begin()+upper+1);

  std::vector<double> y (length, 0.0);


  for (size_t i=0; i<length; ++i)
  {
    const double distance = std::abs(wavelengths[index] - x[i]);

    y[i] = spectrum[lower+i] * aux::normalDistribution(instrument_profile_sigma[index], distance);
  }


  return std::abs(aux::quadratureTrapezoidal(x, y));
}







}
