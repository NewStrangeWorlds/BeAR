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


#include "retrieval.h"


#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>


#include "multinest_parameter.h"
#include "prior.h"
#include "../observations/observations.h"
#include "../forward_model/forward_model.h"
#include "../forward_model/brown_dwarf.h"

#include "../multinest/multinest.h"
#include "../CUDA_kernels/data_management_kernels.h"


namespace helios{



Retrieval::Retrieval(GlobalConfig* global_config) : spectral_grid(global_config)
{

  config = global_config;
  
}




void Retrieval::doRetrieval()
{
  //First, load the observations
  std::string folder = config->retrievalFolder();
  std::string observation_folder = folder;

  //load the observational data
  std::vector<std::string> file_list;
  if ( loadObservationFileList(observation_folder, file_list) == false) return;
  if ( loadObservations(observation_folder, file_list) == false) return;


  std::cout << "Total number of wavelength points: " << spectral_grid.nbSpectralPoints() << "\n";

  
  setPriors();


  //Configure Multinest
  MultinestParameter param;

  param.is = config->multinest_ins;
  param.mmodal = config->multinest_mode_sep;
  param.nlive = config->multinest_nb_living_points;
  param.efr = config->multinest_efficiency;
  param.maxiter = config->multinest_nb_iterations;
  param.resume = config->multinest_resume;
  param.fb = config->multinest_feedback;
  strcpy(param.root, config->retrieval_folder_path.c_str());

  size_t nb_priors = priors.size();

  param.ndims = nb_priors;
  param.nPar = nb_priors;
  param.nClsPar = nb_priors;

  for (size_t i = 0; i < nb_priors; ++i) 
    param.pWrap[i] = 0;
  

  //We give the MultiNest function a pointer to the retrieval class
  //That way, we can access the current Retrieval object in its static member routines
  param.context = this;


  //Initialise the forward model
  double atmos_boundaries [2] = {config->atmosphere_top_pressure, config->atmosphere_bottom_pressure};
  BrownDwarfModel* model = new BrownDwarfModel(this, config->nb_atmosphere_levels, atmos_boundaries);
  
  forward_model = model;



	//Call MultiNest
  if (config->useGPU() == false)
	  nested::run(param.is, param.mmodal, param.ceff, param.nlive, param.tol, param.efr, param.ndims, param.nPar, param.nClsPar,
                param.maxModes, param.updInt, param.Ztol, param.root, param.seed, param.pWrap, param.fb, param.resume,
                param.outfile, param.initMPI, param.logZero, param.maxiter,
                Retrieval::multinestLogLike, Retrieval::multinestDumper,
                param.context);
  else
    nested::run(param.is, param.mmodal, param.ceff, param.nlive, param.tol, param.efr, param.ndims, param.nPar, param.nClsPar,
                param.maxModes, param.updInt, param.Ztol, param.root, param.seed, param.pWrap, param.fb, param.resume,
                param.outfile, param.initMPI, param.logZero, param.maxiter,
                Retrieval::multinestLogLikeGPU, Retrieval::multinestDumper,
                param.context);



  //And finish by deleting stuff
  for (size_t i=0; i<priors.size(); ++i)
    delete priors[i];


  delete forward_model;
}



//load the observational data
//input value is the location of the retrival folder
//return value is a boolean, where false indicates that something went wrong
bool Retrieval::loadObservations(const std::string file_folder, const std::vector<std::string>& file_list)
{
  nb_observations = file_list.size();

  //create the observation object
  observations.resize(nb_observations);

  for (size_t i=0; i<nb_observations; ++i)
  {
    bool init_succes = observations[i].init(this, file_folder + file_list[i]);

    if (!init_succes) return false;
    
    //save all the observations and their errors in a single vector
    observation_data.insert(std::end(observation_data), std::begin(observations[i].flux), std::end(observations[i].flux));
    observation_error.insert(std::end(observation_error), std::begin(observations[i].flux_error), std::end(observations[i].flux_error));
  }


  nb_observation_points = observation_data.size();


  //move the lists to the GPU, if necessary
  if (config->useGPU())
  {
    moveToDevice(observation_data_gpu, observation_data);
    moveToDevice(observation_error_gpu, observation_error);
  }

  

  for (auto & i : observations)
  {
    i.spectral_bands.setLocalIndices();
    i.spectral_bands.setInstrumentProfileFWHW(i.instrument_profile_fwhm);
    i.printObservationDetails();
  }



  //create the vector for the high-res spectrum on the GPU
  if (config->useGPU())
    allocateOnDevice(model_spectrum_gpu, spectral_grid.nbSpectralPoints());


  //if we need to convolve the spectrum with an instrument profile on the GPU, move additional data to the GPU
  //and create the data structures for the convolved spectra of each observation
  if (config->useGPU())
  {
    convolved_spectra.assign(nb_observations, nullptr);
    observation_spectral_indices.assign(nb_observations, nullptr);
    observation_wavelengths.assign(nb_observations, nullptr);
    convolution_start_index.assign(nb_observations, nullptr);
    convolution_end_index.assign(nb_observations, nullptr);
    observation_profile_sigma.assign(nb_observations, nullptr);

   
    for (size_t i=0; i<nb_observations; ++i)
    {
      if (observations[i].instrument_profile_fwhm.size() == 0) continue;
     
      allocateOnDevice(convolved_spectra[i], spectral_grid.nbSpectralPoints());

     
      std::vector<int> band_indices(observations[i].spectral_bands.spectral_indices.begin(), observations[i].spectral_bands.spectral_indices.end());
      moveToDevice(observation_spectral_indices[i], band_indices);

     
      moveToDevice(observation_wavelengths[i], observations[i].spectral_bands.wavelengths);

     
      std::vector<int> start_index(observations[i].spectral_bands.wavelengths.size(), 0);
      std::vector<int> end_index(observations[i].spectral_bands.wavelengths.size(), 0);
     
      for (size_t j=0; j<observations[i].spectral_bands.wavelengths.size(); ++j)
      {
        start_index[j] = observations[i].spectral_bands.convolution_quadrature_intervals[j][0];
        end_index[j] = observations[i].spectral_bands.convolution_quadrature_intervals[j][1];
      }
      
     
      moveToDevice(convolution_start_index[i], start_index);
      moveToDevice(convolution_end_index[i], end_index);

      moveToDevice(observation_profile_sigma[i], observations[i].spectral_bands.instrument_profile_sigma);
    }

  }



  //for the GPU we move all sub-band indices into one vector and
  //additionally store also the number of points per sub-band
  if (config->useGPU())
  {
    std::vector<int> band_indices_host;
    std::vector<int> band_sizes_host;
    std::vector<int> band_start_index_host;

    for (size_t i=0; i<nb_observations; ++i)
    {

      for (size_t j=0; j<observations[i].spectral_bands.nbBands(); ++j)
      {

        band_start_index_host.push_back(band_indices_host.size());

        band_indices_host.insert(std::end(band_indices_host),
                                 std::begin(observations[i].spectral_bands.band_spectral_indices[j]),
                                 std::end(observations[i].spectral_bands.band_spectral_indices[j]));

        band_sizes_host.push_back(observations[i].spectral_bands.band_spectral_indices[j].size());

        
        //assign a high-res spectrum to each band
        //if a band from a certain observation requires convolution, assign the pointer to the convolved spectrum
        if (observations[i].instrument_profile_fwhm.size() == 0)
          band_spectrum_id.push_back(model_spectrum_gpu);
        else
          band_spectrum_id.push_back(convolved_spectra[i]);


        nb_total_bands++;
      }

    }



    moveToDevice(band_indices_gpu, band_indices_host);
    moveToDevice(band_sizes_gpu, band_sizes_host);
    moveToDevice(band_start_index_gpu, band_start_index_host);

  }


  std::cout << "loaded observations: \n";

  for (size_t i=0; i<nb_observations; ++i)
    std::cout << "name " << observations[i].observationName() << "\n";

  return true;
}




//load the observational file list
//input value is the location of the retrival folder
//return value is a boolean, where false indicates that something went wrong
bool Retrieval::loadObservationFileList(const std::string file_folder, std::vector<std::string>& file_list)
{

  //we are look for the file observations.list in the folder
  std::string file_name = file_folder + "observations.list";
  std::fstream file;

  file.open(file_name.c_str(), std::ios::in);

  //if we don't find it, we return false
  if (file.fail())
  {
    std::cout << "Could not open observation list file: " << file_name << "\n";
    return false;
  }


  //read the list of observation data files
  std::string observation_file;

  while (file >> observation_file)
    file_list.push_back(observation_file);

  nb_observations = file_list.size();

  return true;
}



//the destructor
//deletes all remaining data on the GPU
Retrieval::~Retrieval()
{

  if (config->useGPU())
  {
    deleteFromDevice(observation_data_gpu);
    deleteFromDevice(observation_error_gpu);

    deleteFromDevice(band_sizes_gpu);
    deleteFromDevice(band_indices_gpu);
    deleteFromDevice(band_start_index_gpu);
  }


}


}
