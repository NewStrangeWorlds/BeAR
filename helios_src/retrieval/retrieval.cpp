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



namespace helios{


bool stop_model = false;

void signalHandler(int sig) 
{
  std::cout << "Received signal " << sig << "\n"; 
  
  if (sig == SIGCONT) 
    stop_model = true;

}


Retrieval::Retrieval(GlobalConfig* global_config) 
  : spectral_grid(global_config)
{

  config = global_config;

  std::signal(SIGCONT, signalHandler);

}




bool Retrieval::doRetrieval()
{
  std::string folder = config->retrieval_folder_path;
  std::string observation_folder = folder;


  //try to initialise the model
  //if there is an error, we exit the retrieval
  try
  {
    std::vector<std::string> file_list;
    loadObservationFileList(observation_folder, file_list);
    loadObservations(observation_folder, file_list);

    std::cout << "\nTotal number of wavelength points: " << spectral_grid.nbSpectralPoints() << "\n\n";

    forward_model = selectForwardModel(config->forward_model_type);
  }
  catch(std::runtime_error& e) 
  {
    std::cout << e.what() << std::endl;
    return false;
  }


  setAdditionalPriors();

  priors.printInfo();


  //Configure Multinest
  MultinestParameter param(config);

  size_t nb_priors = priors.number();

  param.ndims = nb_priors;
  param.nPar = nb_priors;
  param.nClsPar = nb_priors;

  for (size_t i = 0; i < nb_priors; ++i) 
    param.pWrap[i] = 0;

  //We give the MultiNest function a pointer to the retrieval class
  //That way, we can access the current retrieval object in its static member routines
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
    std::vector<double>::iterator it = std::min_element(std::begin(observation_error), std::end(observation_error));
    const double error_min = std::log10(0.1 * *it * *it);

    it = std::max_element(std::begin(observation_error), std::end(observation_error));
    const double error_max = std::log10(100.0 * *it * *it);

    priors.add(
      std::vector<std::string>{std::string("uniform")}, 
      std::vector<std::string>{std::string("error exponent")}, 
      std::vector<std::vector<double>>{std::vector<double> {error_min, error_max}});
  }
}



Retrieval::~Retrieval()
{
  if (config->use_gpu)
  {
    deleteFromDevice(model_spectrum_gpu);

    deleteFromDevice(observation_data_gpu);
    deleteFromDevice(observation_error_gpu);
    deleteFromDevice(observation_likelihood_weight_gpu);
  }
}


}
