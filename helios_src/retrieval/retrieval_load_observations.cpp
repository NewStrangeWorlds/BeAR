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


#include "retrieval.h"


#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>


#include "../observations/observations.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../additional/exceptions.h"



namespace helios{


//load the observational data
//input value is the location of the retrieval folder
void Retrieval::loadObservations(const std::string file_folder, const std::vector<std::string>& file_list)
{
  nb_observations = file_list.size();

  //create the observation object
  //observations.resize(nb_observations);

  observations.assign(nb_observations, Observation(config, &spectral_grid));

  for (size_t i=0; i<nb_observations; ++i)
  {
    observations[i].init(file_folder + file_list[i]);
    
    //save all the observations and their errors in a single vector
    observation_data.insert(std::end(observation_data), std::begin(observations[i].flux), std::end(observations[i].flux));
    observation_error.insert(std::end(observation_error), std::begin(observations[i].flux_error), std::end(observations[i].flux_error));
    observation_likelihood_weight.insert(std::end(observation_likelihood_weight), std::begin(observations[i].likelihood_weight), std::end(observations[i].likelihood_weight));
  }


  nb_observation_points = observation_data.size();

  for (auto & i : observations)
    nb_total_bands += i.flux.size();


  //move the lists to the GPU, if necessary
  if (config->use_gpu)
  {
    moveToDevice(observation_data_gpu, observation_data);
    moveToDevice(observation_error_gpu, observation_error);
    moveToDevice(observation_likelihood_weight_gpu, observation_likelihood_weight);
  }


  for (auto & i : observations)
  {
    i.spectral_bands.setLocalIndices();
    i.spectral_bands.setInstrumentProfileFWHW(i.instrument_profile_fwhm);
    i.setFilterResponseFunction();
    i.printObservationDetails();
    i.spectral_bands.initDeviceMemory();
    i.initDeviceMemory();
  }


  //create the vector for the high-res spectrum on the GPU
  if (config->use_gpu)
    allocateOnDevice(model_spectrum_gpu, spectral_grid.nbSpectralPoints());

  std::cout << "loaded observations: \n";

  for (size_t i=0; i<nb_observations; ++i)
    std::cout << "name " << observations[i].observationName() << "\n";
}




//load the observational file list
//input value is the location of the retrival folder
void Retrieval::loadObservationFileList(const std::string file_folder, std::vector<std::string>& file_list)
{
  //we are look for the file observations.list in the folder
  std::string file_name = file_folder + "observations.list";
  std::fstream file;

  file.open(file_name.c_str(), std::ios::in);

  
  if (file.fail())
    throw ExceptionFileNotFound(std::string ("Retrieval::loadObservationFileList"), file_name);


  //read the list of observation data files
  std::string observation_file;

  while (file >> observation_file)
    file_list.push_back(observation_file);

  nb_observations = file_list.size();


  if (nb_observations == 0)
  {
    std::string error_message = "No observations found in observations.list file.\n";
    throw ExceptionInvalidInput(std::string ("Retrieval::loadObservationFileList"), error_message);
  }
}



}
