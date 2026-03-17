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
#include <algorithm>

#include "retrieval.h"

#include "../observations/observations.h"
#include "../observations/highres_observation.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../additional/exceptions.h"


namespace bear{



void Retrieval::setObservations(
  const std::vector<ObservationInput>& observation_input)
{
  nb_observations = observation_input.size();

  observations.assign(nb_observations, Observation(config, &spectral_grid));

  for (size_t i=0; i<nb_observations; ++i)
    observations[i].init(observation_input[i]);

  spectral_grid.sampleSpectralGrid(observations);

  for (auto & i : observations)
  {
    i.spectral_bands.setInstrumentProfileFWHW(i.instrument_profile_fwhm);
    i.setFilterResponseFunction();
    i.printObservationDetails();
    i.spectral_bands.initDeviceMemory();
    i.initDeviceMemory();
  }

  std::cout << "loaded observations: \n";

  for (size_t i=0; i<nb_observations; ++i)
    std::cout << "name " << observations[i].observationName() << "\n";
}



void Retrieval::loadObservations(
  const std::string file_folder, 
  const std::vector<std::string>& file_list,
  const std::vector<std::string>& modifier_list)
{
  nb_observations = file_list.size();

  observations.assign(nb_observations, Observation(config, &spectral_grid));

  for (size_t i=0; i<nb_observations; ++i)
    observations[i].init(file_folder + file_list[i], modifier_list[i]);
    

  spectral_grid.sampleSpectralGrid(observations);


  for (auto & i : observations)
  {
    i.spectral_bands.setInstrumentProfileFWHW(i.instrument_profile_fwhm);
    i.setFilterResponseFunction();
    i.printObservationDetails();
    i.spectral_bands.initDeviceMemory();
    i.initDeviceMemory();
  }

  std::cout << "loaded observations: \n";

  for (size_t i=0; i<nb_observations; ++i)
    std::cout << "name " << observations[i].observationName() << "\n";
}




//load the observational file list
//input value is the location of the retrival folder
void Retrieval::loadObservationFileList(
  const std::string file_folder,
  std::vector<std::string>& file_list,
  std::vector<std::string>& modifier_list)
{
  //we are look for the file observations.list in the folder
  std::string file_name = file_folder + "observations.list";
  std::fstream file;

  file.open(file_name.c_str(), std::ios::in);


  if (file.fail())
    throw FileNotFound(std::string ("Retrieval::loadObservationFileList"), file_name);

  std::string line;

  while (std::getline(file, line))
  {
    if (line.empty() || line[0] == '#') continue;

    std::string observation_file = "";
    std::string observation_modifier = "";

    std::stringstream ss(line);

    ss >> observation_file >> observation_modifier;

    if (observation_modifier == "high_resolution")
    {
      highres_file_list.push_back(observation_file);
    }
    else
    {
      file_list.push_back(observation_file);
      modifier_list.push_back(observation_modifier);
    }
  }

  nb_observations = file_list.size();

  if (nb_observations == 0 && highres_file_list.empty())
  {
    std::string error_message = "No observations found in observations.list file.\n";
    throw InvalidInput(std::string ("Retrieval::loadObservationFileList"), error_message);
  }
}


void Retrieval::loadHighResObservations(const std::string& file_folder)
{
  if (highres_file_list.empty()) return;

  nb_highres_observations = highres_file_list.size();
  has_highres_observations = true;

  highres_observations.resize(nb_highres_observations);

  for (size_t i = 0; i < nb_highres_observations; ++i)
    highres_observations[i].init(file_folder + highres_file_list[i]);

  // Determine overall wavelength range from all high-res observations
  double wl_min = 1e30;
  double wl_max = 0;

  for (const auto& obs : highres_observations)
  {
    wl_min = std::min(wl_min, obs.wavelengthMin());
    wl_max = std::max(wl_max, obs.wavelengthMax());
  }

  // Create the high-res spectral grid using spectral_resolution_highres.
  // Extend the wavelength range by a velocity margin to accommodate
  // Doppler shifts up to ~300 km/s without edge pixels falling outside
  // the model range and returning zero from interpolation.
  const double doppler_margin = 300.0 / 299792.458;  // ~300 km/s in v/c

  double wl_min_ext = wl_min * (1.0 - doppler_margin);
  double wl_max_ext = wl_max * (1.0 + doppler_margin);

  // Temporarily override config to use the high-res resolution
  double saved_resolution = config->spectral_resolution;
  unsigned int saved_disc = config->spectral_disecretisation;

  config->spectral_resolution = config->spectral_resolution_highres;
  config->spectral_disecretisation = 2;  // constant resolving power mode

  // Convert nm to microns (SpectralGrid works in microns internally)
  spectral_grid_highres = std::make_unique<SpectralGrid>(
    config, wl_min_ext * 1e-3, wl_max_ext * 1e-3);

  // Restore original config values
  config->spectral_resolution = saved_resolution;
  config->spectral_disecretisation = saved_disc;

  std::cout << "High-res spectral grid: "
            << spectral_grid_highres->nbSpectralPoints()
            << " points, R = " << config->spectral_resolution_highres
            << ", range: " << wl_min << " - " << wl_max << " nm\n";

  // GPU memory for high-res observations is initialized later,
  // after the likelihood mode has been determined (see retrieval init).
}



}
