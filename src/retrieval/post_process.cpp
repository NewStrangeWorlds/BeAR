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
#include <cstdio>

#include "post_process.h"

#include "../observations/observations.h"
#include "../forward_model/forward_model.h"
#include "../additional/physical_const.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../additional/exceptions.h"
#include "retrieval.h"


namespace bear{


PostProcess::PostProcess(GlobalConfig* global_config) 
  : Retrieval(global_config, std::string("postprocess_spectrum_data.dat"))
{
  
  config = global_config;
  
}



bool PostProcess::run()
{
  std::string folder = config->retrieval_folder_path;
  std::string observation_folder = folder;
  bool delete_sampler_files = false;

  try
  {
    readPosteriorData();
    
    forward_model->postProcess(model_parameter, best_fit_model, delete_sampler_files);
  }
  catch(std::runtime_error& e) 
  {
    std::cout << e.what() << std::endl;
   
    return false;
  }
  
  if (delete_sampler_files)
  {
    MultinestParameter param(config);
    deleteSamplerFiles(param.unused_posterior_files);
  }

  return true;
}



void PostProcess::deleteSamplerFiles(const std::vector<std::string>& file_list)
{ 
  
  for (auto & f : file_list)
  {
    std::string file_name = config->retrieval_folder_path + f;
    std::remove(file_name.c_str());
  }

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

  //scale the model parameters with their units
  for (auto & m : model_parameter)
    for (size_t i=0; i<priors.distributions.size(); ++i)
    {
      m[i] = priors.distributions[i]->applyParameterUnit(m[i]);
    }

  //find best-fit model
  best_fit_model = 0;
  best_log_like = log_like[0];

  for (size_t i=1; i<log_like.size(); ++i)
    if (log_like[i] > best_log_like)
    {
      best_fit_model = i;
      best_log_like = log_like[i];
    }

  std::cout << "\nBest-fit model: " << best_fit_model << "\t ln(Z): " << best_log_like << "\n";
}



PostProcess::~PostProcess()
{
  
}


}
