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


#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <iomanip>
#include <sstream>

#include "star_file_spectrum.h"

#include "../../spectral_grid/spectral_grid.h"
#include "../../additional/physical_const.h"
#include "../../additional/aux_functions.h"
#include "../../additional/exceptions.h"


namespace bear{


StarSpectrumFile::StarSpectrumFile(
  const std::string file_path,
  SpectralGrid* spectral_grid_)
  : spectral_grid(spectral_grid_)
{
  nb_parameters = 0;
  
  std::vector<double> spectrum_file;
  std::vector<double> wavelength_file;

  readSpectrum(file_path, spectrum_file, wavelength_file);

  spectrum = spectral_grid->interpolateToWavelengthGrid(wavelength_file, spectrum_file, false);
}


StarSpectrumFile::StarSpectrumFile(
  const std::vector<double>& wavelengths,
  const std::vector<double>& flux_,
  SpectralGrid* spectral_grid_)
  : spectral_grid(spectral_grid_)
{
  nb_parameters = 0;
  
  std::vector<double> flux = flux_;

  //convert from W m-2 mu-1 to W m-2 cm
  for (size_t i=0; i<flux.size(); ++i)
    flux[i] = flux[i]*wavelengths[i]*wavelengths[i]/10000.;

  spectrum = spectral_grid->interpolateToWavelengthGrid(wavelengths, flux, false);
}


std::vector<double> StarSpectrumFile::calcFlux(
  const std::vector<double>& parameter)
{

  return spectrum;

}



void StarSpectrumFile::readSpectrum(
  const std::string file_path,
  std::vector<double>& spectrum_file,
  std::vector<double>& wavelength_file)
{
  std::fstream file;

  file.open(file_path.c_str(), std::ios::in);


  if (file.fail())
    throw FileNotFound(std::string ("StellarSpectrum::readSpectrum"), file_path);


  std::cout << "Reading stellar spectrum file " << file_path << "\n";
  
  wavelength_file.reserve(5000000);
  spectrum_file.reserve(5000000);

  std::string line;

  while (std::getline(file, line))
  {
    std::stringstream line_stream(line);

    double wavelength_in;
    double spectrum_in;

    if (!(line_stream >> wavelength_in >> spectrum_in)) continue;

    wavelength_file.push_back(wavelength_in);  
    spectrum_file.push_back(spectrum_in);
  }

  file.close();

  wavelength_file.shrink_to_fit(); 
  spectrum_file.shrink_to_fit();

  //convert from W m-2 mu-1 to W m-2 cm
  for (size_t i=0; i<spectrum_file.size(); ++i)
    spectrum_file[i] = spectrum_file[i]*wavelength_file[i]*wavelength_file[i]/10000.;
}


}