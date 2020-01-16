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


#include "post_process.h"


#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>


#include "../observations/observations.h"
#include "../forward_model/forward_model.h"
#include "../forward_model/brown_dwarf.h"
#include "../additional/physical_const.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../CUDA_kernels/band_integration_kernels.h"
#include "../forward_model/piecewise_poly.h"
#include "retrieval.h"



namespace helios{



PostProcess::PostProcess(GlobalConfig* global_config) : Retrieval(global_config)
{
  
  config = global_config;
  
}



void PostProcess::doRetrieval()
{

  std::string folder = config->retrievalFolder();
  std::string observation_folder = folder;

  std::vector<std::string> file_list;
  if ( loadObservationFileList(observation_folder, file_list) == false) return;
  
  //if we do postprocessing, we need to read in the file that describes the maximum wavelength range
  //spectra will be generated for
  //this is necessary to obtain an estimate for the effective temperature 
  file_list.push_back("postprocess_spectrum_data.dat");
  
  if ( loadObservations(observation_folder, file_list) == false) return;

  std::cout << "\nTotal number of wavelength points: " << spectral_grid.nbSpectralPoints() << "\n";


  
  //Initialise the forward model
  double atmos_boundaries [2] = {100, 1e-3};
  BrownDwarfModel* model = new BrownDwarfModel(this, config->nb_atmosphere_levels, atmos_boundaries);
  
  forward_model = model;

  readPosteriorData();


  std::vector< std::vector<double> > model_spectrum_bands;  

  postProcessSpectra(model_spectrum_bands);

  std::vector<double> effective_temperatures;

  postProcessEffectiveTemperatures(model_spectrum_bands, effective_temperatures);


  saveOutput(model_spectrum_bands, effective_temperatures);

  
  postProcessTemperatures();


  delete model;
}



//Read the posterior data from the Multinest output file
bool PostProcess::readPosteriorData()
{ 
  //first, we load the results from the written file
  std::fstream file;

  std::string file_name = config->retrieval_folder_path.c_str();
  file_name += "post_equal_weights.dat";

  file.open(file_name.c_str(),std::ios::in);


  if (file.fail())
  {
    std::cout << "Could not open posterior file: " << file_name << "\n";
    return false;
  }


  std::string line;
  size_t nb_param = 0; 

  //first ne need to know the number of parameters
  std::getline(file, line);

  std::stringstream line_stream(line);
  
  double val = 0;
  while (line_stream >> val)
    nb_param++;

  //the last column is the likelihood
  nb_param--;


  //return to the beginning of the file     
  file.seekg(0, std::ios::beg);


  //read the posteriors and likelihoods
  while(std::getline(file, line))
  {
    if (file.eof()) break;
    if (line.size() < 2) break;

    std::stringstream line_stream(line);

    std::vector<double> param(nb_param, 0.0);
    double log_z;

    for (size_t i=0; i<nb_param; ++i)
      line_stream >> param[i];

    line_stream >> log_z;


    model_parameter.push_back(param);
    log_like.push_back(log_z);
  }

  file.close();


  //find best-fit model
  best_fit_model = 0;
  best_log_like = log_like[0];

  for (size_t i=1; i<log_like.size(); ++i)
    if (log_like[i] > best_log_like)
    {

      best_fit_model = i;
      best_log_like = log_like[i];

    }

  std::cout << "best fit model " << best_fit_model << "\t" << best_log_like << "\n";


  return true;
}



//Save the post-process output to the disk
void PostProcess::saveOutput(const std::vector< std::vector<double> >& model_spectrum_bands, const std::vector<double>& effective_temperatures)
{

  const size_t nb_models = model_parameter.size();

  //write the spectra to files
  unsigned int band_index = 0;

  for (size_t j=0; j<nb_observations; ++j)
  {
    std::string file_name = config->retrievalFolder() + "/spectrum_lr_o" + std::to_string(j) + ".dat";

    std::fstream file(file_name.c_str(), std::ios::out);

    for (size_t i=0; i<observations[j].spectral_bands.nbBands(); ++i)
    {

      file << std::setprecision(10) << std::scientific << observations[j].spectral_bands.band_centers_wavelength[i];
      
      for (size_t k=0; k<nb_models; ++k)
        file << "\t" << model_spectrum_bands[k][band_index + i];
     
      file << "\n";
    }


    band_index += observations[j].spectral_bands.nbBands();
  }


  //save the effective temperatures
  std::string file_name = config->retrievalFolder() + "/effective_temperatures.dat";
  
  std::fstream file(file_name.c_str(), std::ios::out);

  for (size_t k=0; k<nb_models; ++k)
    file << std::setprecision(10) << std::scientific << effective_temperatures[k] << "\n";
 

}





PostProcess::~PostProcess()
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
