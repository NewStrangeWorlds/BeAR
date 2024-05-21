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


#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "observations.h"

#include "../spectral_grid/spectral_grid.h"
#include "../additional/quadrature.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"
#include "../additional/aux_functions.h"
#include "../CUDA_kernels/data_management_kernels.h"


namespace helios{


std::vector<std::vector<double>> Observation::readFilterResponseFunction(
  const std::string& file_path)
{
  std::fstream file;
  
  file.open(file_path.c_str(), std::ios::in);


  if (file.fail())
    throw FileNotFound(std::string ("SpectralBands::readFilterResponseFunction"), file_path);
 
  
  std::cout << "\nReading filter response file " << file_path << "\n";
  
  std::vector<double> wavelength;  wavelength.reserve(10000);
  std::vector<double> response_function;  response_function.reserve(10000);

  std::string line; 

  //header
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);

  std::stringstream line_stream(line);
  std::string detector_type;

  line_stream >> detector_type;
  std::cout << "Detector type: " << detector_type << "\n";

  if (detector_type == "Energy")
    detector_type = "energy";

  if (detector_type == "Photon")
    detector_type = "photon";

  if (detector_type != "energy" && detector_type != "photon")
  {
    std::string error_message = "Unsupported detector type *" + detector_type + "* in file: " + file_path + "\n";
    throw InvalidInput(std::string ("Observation::readFilterResponseFunction"), error_message);
  }

  filter_detector_type = detector_type;

  //read in the response function
  while (std::getline(file, line))
  {
    std::stringstream line_stream(line);

    double wavelength_in;
    double spectrum_in;

    if (!(line_stream >> wavelength_in >> spectrum_in)) continue;

    wavelength.push_back(wavelength_in);  
    response_function.push_back(spectrum_in);
  }

  file.close();

  wavelength.shrink_to_fit(); response_function.shrink_to_fit();


  if (wavelength.front() < wavelength.back())
  {
    std::reverse(wavelength.begin(), wavelength.end());
    std::reverse(response_function.begin(), response_function.end());
  }


  std::vector<std::vector<double>> filter_response_data = {wavelength, response_function};


  return filter_response_data;
}




void Observation::setFilterResponseFunction()
{
  if (filter_response_file.size() == 0) return;

  filter_response = spectral_grid->interpolateToWavelengthGrid(filter_response_file[0], filter_response_file[1], false);

  filter_response_weight.assign(filter_response.size(), 1.0);


  if (filter_detector_type == "photon")
    for (size_t i=0; i<spectral_grid->nbSpectralPoints(); ++i)
      filter_response_weight[i] = spectral_grid->wavelength_list[i];


  if (config->use_gpu)
  {
    moveToDevice(filter_response_gpu, filter_response);
    moveToDevice(filter_response_weight_gpu, filter_response_weight);
  }
  
  filter_response_normalisation = filterResponseNormalisation(filter_response_file[0], filter_response_file[1]);
}



double Observation::filterResponseNormalisation(const std::vector<double>& wavelength, const std::vector<double>& response_function)
{
  std::vector<double> y = response_function;
  
  if (filter_detector_type == "photon")
    for (size_t i=0; i<y.size(); ++i)
      y[i] *= wavelength[i];

  double normalisation = aux::quadratureTrapezoidal(wavelength, y);

  if (normalisation < 0) normalisation = std::abs(normalisation);

  return normalisation;
}




std::vector<double> Observation::applyFilterResponseFunction(const std::vector<double>& spectrum)
{
  if (filter_response.size() == 0) return spectrum;


  std::vector<double> filter_spectrum(spectrum.size(), 0.0);

  for (size_t i=0; i<spectral_grid->nbSpectralPoints(); ++i)
  {

    filter_spectrum[i] = spectrum[i] * filter_response[i] * filter_response_weight[i] / filter_response_normalisation;

  }


  return filter_spectrum;
}



}
