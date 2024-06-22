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
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>

#include "post_process.h"

#include "../observations/observations.h"
#include "../forward_model/forward_model.h"
#include "../additional/physical_const.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../additional/exceptions.h"
#include "retrieval.h"


namespace bear{



PostProcess::PostProcess(GlobalConfig* global_config) : Retrieval(global_config)
{
  
  config = global_config;
  
}



bool PostProcess::doRetrieval()
{
  std::string folder = config->retrieval_folder_path;
  std::string observation_folder = folder;


  //try to initialise the model
  //if there is an error, we exit the post process
  try
  {
    std::vector<std::string> file_list, modifier_list;
    loadObservationFileList(
      observation_folder,
      file_list,
      modifier_list);
  
    //if we do postprocessing, we may need to read in the file that describes the maximum wavelength range
    //spectra will be generated for
    //this is necessary to obtain an estimate for the effective temperature
    std::string postprocess_spectrum_data = config->retrieval_folder_path + "postprocess_spectrum_data.dat";
    std::fstream file(postprocess_spectrum_data.c_str(), std::ios::in);

    if (!file.fail())
    {
      file_list.push_back("postprocess_spectrum_data.dat");
      modifier_list.push_back("none");
      file.close(); 
    }

    loadObservations(
      observation_folder,
      file_list,
      modifier_list);

    std::cout << "\nTotal number of wavelength points: " << spectral_grid.nbSpectralPoints() << "\n\n";

    //Initialise the forward model
    forward_model = selectForwardModel(config->forward_model_type);

    readPosteriorData();
  }
  catch(std::runtime_error& e) 
  {
    std::cout << e.what() << std::endl;
    return false;
  } 


  std::vector< std::vector<double> > model_spectrum_bands;  

  postProcessSpectra(model_spectrum_bands);

  saveOutput(model_spectrum_bands);


  //call the post postprocess method of the forward model
  forward_model->postProcess(model_parameter, model_spectrum_bands, best_fit_model);


  delete forward_model;


  return true;
}


void PostProcess::readPosteriorData()
{ 
  //first, we load the results from the written file
  std::fstream file;

  std::string file_name = config->retrieval_folder_path.c_str();
  file_name += "post_equal_weights.dat";

  file.open(file_name.c_str(),std::ios::in);


  if (file.fail())
    throw FileNotFound(std::string ("PostProcess::readPosteriorData"), file_name);


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
    if (line.size() < 5) break;

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
}


void PostProcess::saveOutput(const std::vector< std::vector<double> >& model_spectrum_bands)
{
  const size_t nb_models = model_parameter.size();

  //write the spectra to files
  unsigned int band_index = 0;

  for (size_t j=0; j<nb_observations; ++j)
  {
    std::string observation_name = observations[j].observationName();
    std::replace(observation_name.begin(), observation_name.end(), ' ', '_'); 
    
    std::string file_name = config->retrieval_folder_path + "/spectrum_post_" + observation_name + ".dat"; 
  

    std::fstream file(file_name.c_str(), std::ios::out);

    for (size_t i=0; i<observations[j].spectral_bands.nbBands(); ++i)
    {

      file << std::setprecision(10) << std::scientific << observations[j].spectral_bands.center_wavelengths[i];
      
      for (size_t k=0; k<nb_models; ++k)
        file << "\t" << model_spectrum_bands[k][band_index + i];
     
      file << "\n";
    }

    file.close();

    band_index += observations[j].spectral_bands.nbBands();
  }

}



PostProcess::~PostProcess()
{
  if (config->use_gpu)
  {
    deleteFromDevice(observation_data_gpu);
    deleteFromDevice(observation_error_gpu);
  }

}


}
