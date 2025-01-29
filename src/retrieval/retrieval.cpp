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


#include "retrieval.h"


#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <csignal>
#include <cstdlib>


#include "multinest_parameter.h"
#include "priors.h"
#include "../observations/observations.h"
#include "../forward_model/forward_model.h"

#include "../../_deps/multinest-src/MultiNest_v3.12_CMake/multinest/include/multinest.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../additional/exceptions.h"



namespace bear{


bool stop_model = false;

void signalHandler(int sig) 
{
  std::cout << "Received signal " << sig << "\n"; 
  
  if (sig == SIGCONT) 
    stop_model = true;

}


Retrieval::Retrieval(GlobalConfig* global_config) 
  : Retrieval(global_config, std::string(""))
{
  config = global_config;

  std::signal(SIGCONT, signalHandler);
}



Retrieval::Retrieval(
  GlobalConfig* global_config,
  const std::string additional_observation_file) 
  : spectral_grid(global_config)
{
  config = global_config;

  std::signal(SIGCONT, signalHandler);

  std::string folder = config->retrieval_folder_path;
  std::string observation_folder = folder;

  //try to initialise the model
  //if there is an error, we exit the retrieval
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
    if (additional_observation_file.size() > 0)
    {
      std::string postprocess_spectrum_data = config->retrieval_folder_path + additional_observation_file;
      std::fstream file(postprocess_spectrum_data.c_str(), std::ios::in);

      if (!file.fail())
      {
        file_list.push_back(additional_observation_file);
        modifier_list.push_back("none");
        file.close(); 
      }
    }

    loadObservations(
      observation_folder,
      file_list,
      modifier_list);
  
    std::cout << "\nTotal number of wavelength points: " << spectral_grid.nbSpectralPoints() << "\n\n";

    forward_model = selectForwardModel(config->forward_model_type);
  }
  catch(std::runtime_error& e) 
  {
    std::cout << e.what() << std::endl;
    exit(1);
  }

  setAdditionalPriors();

  priors.printInfo();
}




bool Retrieval::run()
{
  //Configure Multinest
  MultinestParameter param(config);

  size_t nb_priors = priors.number();

  param.ndims = nb_priors;
  param.nPar = nb_priors;
  param.nClsPar = nb_priors;

  for (size_t i = 0; i < nb_priors; ++i) 
    param.pWrap[i] = 0;

  //We give the MultiNest function a pointer to the retrieval class
  //That way, we can access the current retrieval object and its static member routines
  param.context = this;


  //Call MultiNest
  if (config->use_gpu == false)
    nested::run(
      param.is,
      param.mmodal,
      param.ceff,
      param.nlive,
      param.tol,
      param.efr,
      param.ndims,
      param.nPar,
      param.nClsPar,
      param.maxModes,
      param.updInt,
      param.Ztol,
      param.root,
      param.seed,
      param.pWrap,
      param.fb,
      param.resume,
      param.outfile,
      param.initMPI,
      param.logZero,
      param.maxiter,
      Retrieval::multinestLogLike,
      Retrieval::multinestDumper,
      param.context);
  else
    nested::run(
      param.is,
      param.mmodal,
      param.ceff,
      param.nlive,
      param.tol,
      param.efr,
      param.ndims,
      param.nPar,
      param.nClsPar,
      param.maxModes,
      param.updInt,
      param.Ztol,
      param.root,
      param.seed,
      param.pWrap,
      param.fb,
      param.resume,
      param.outfile,
      param.initMPI,
      param.logZero,
      param.maxiter,
      Retrieval::multinestLogLikeGPU,
      Retrieval::multinestDumper,
      param.context);


  delete forward_model;


  return true;
}



void Retrieval::setAdditionalPriors()
{
  if (config->use_error_inflation)
  { 
    //this creates the prior distribution for the error exponent
    //first, we need to find the minimum and maximum values of the observational data errors
    double error_max = 0;

    for (auto & obs : observations)
    {
      double obs_error_max = *std::max_element(
        std::begin(obs.data_error), 
        std::end(obs.data_error));

      if (obs_error_max > error_max)
        error_max = obs_error_max;
    }
 
    double error_min = error_max;
    
    for (auto & obs : observations)
    {
      double obs_error_min = *std::min_element(
        std::begin(obs.data_error), 
        std::end(obs.data_error));
      
      if (obs_error_min < error_min)
        error_min = obs_error_min;
    }

    error_min = std::log10(0.1 * error_min * error_min);
    error_max = std::log10(100.0 * error_max * error_max);

    priors.add(
      std::vector<std::string>{std::string("uniform")}, 
      std::vector<std::string>{std::string("error exponent")}, 
      std::vector<std::vector<std::string>>{std::vector<std::string> {std::to_string(error_min), std::to_string(error_max)}});
  }
}



Retrieval::~Retrieval()
{
  
}


}
