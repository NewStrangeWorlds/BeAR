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


#ifndef _stellar_spectrum_h
#define _stellar_spectrum_h

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>


#include "../spectral_grid/spectral_grid.h"
#include "../observations/observations.h"
#include "../additional/exceptions.h"
#include "../CUDA_kernels/data_management_kernels.h"


namespace helios {


class StellarSpectrum {
  public:
    StellarSpectrum (
      const std::string file_path,
      const bool use_gpu,
      SpectralGrid* spectral_grid,
      std::vector<Observation>& observations);
    ~StellarSpectrum();

    std::vector<double> flux;
    double* flux_dev = nullptr;
  private:
    void readSpectrum(
      const std::string file_path,
      std::vector<double>& spectrum_file,
      std::vector<double>& wavelength_file);
    void binSpectrum(
      std::vector<double>& spectrum_file,
      std::vector<double>& wavelength_file,
      SpectralGrid* spectral_grid,
      std::vector<Observation>& observations);
};



inline StellarSpectrum::StellarSpectrum (
  const std::string file_path,
  const bool use_gpu,
  SpectralGrid* spectral_grid,
  std::vector<Observation>& observations)
{
  std::vector<double> spectrum_file;
  std::vector<double> wavelength_file;

  readSpectrum(file_path, spectrum_file, wavelength_file);

  binSpectrum(spectrum_file, wavelength_file, spectral_grid, observations);

  if (use_gpu)
    moveToDevice(flux_dev, flux);
}



inline StellarSpectrum::~StellarSpectrum()
{
  if (flux_dev != nullptr)
    deleteFromDevice(flux_dev);
}



inline void StellarSpectrum::readSpectrum(
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


inline void StellarSpectrum::binSpectrum(
  std::vector<double>& spectrum_file,
  std::vector<double>& wavelength_file,
  SpectralGrid* spectral_grid,
  std::vector<Observation>& observations)
{ 
  std::vector<double> flux_high_res = 
    spectral_grid->interpolateToWavelengthGrid(wavelength_file, spectrum_file, false);
  
  size_t nb_observation_points = 0;

  for (auto & i : observations)
    nb_observation_points += i.nbPoints();

  flux.assign(nb_observation_points, 0.0);

  std::vector<double>::iterator it = flux.begin();

  for (size_t i=0; i<observations.size(); ++i)
  {
    const bool is_flux = true;
    std::vector<double> flux_bands = observations[i].processModelSpectrum(flux_high_res, is_flux);

    std::copy(flux_bands.begin(), flux_bands.end(), it);
    it += flux_bands.size();
  }
  std::cout << "obs " << nb_observation_points << "\n";
  for (size_t i=0; i<flux.size(); ++i)
    std::cout << i << "\t" << flux[i] << "\n";
}


}


#endif

