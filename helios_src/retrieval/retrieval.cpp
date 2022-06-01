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
#include <csignal>
#include <cstdlib>


#include "multinest_parameter.h"
#include "prior.h"
#include "../observations/observations.h"
#include "../forward_model/forward_model.h"
#include "../forward_model/brown_dwarf/brown_dwarf.h"
#include "../forward_model/secondary_eclipse/secondary_eclipse.h"

#include "../../multinest/multinest.h"
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


Retrieval::Retrieval(GlobalConfig* global_config) : spectral_grid(global_config)
{

  config = global_config;

  std::signal(SIGCONT, signalHandler);

}




bool Retrieval::doRetrieval()
{
  std::string folder = config->retrieval_folder_path;
  std::string observation_folder = folder;


  BrownDwarfModel* model = nullptr;
  //SecondaryEclipseModel* model = nullptr;


  //try to initialise the model
  //if there is an error, we exit the retrieval
  try
  {
    //load the observational data
    std::vector<std::string> file_list;
    loadObservationFileList(observation_folder, file_list);
    loadObservations(observation_folder, file_list);


    std::cout << "\nTotal number of wavelength points: " << spectral_grid.nbSpectralPoints() << "\n\n";
  

    //Initialise the forward model
    model = new BrownDwarfModel(this, BrownDwarfConfig (config->retrieval_folder_path));
    //model = new SecondaryEclipseModel(this, SecondaryEclipseConfig (config->retrieval_folder_path));
  }
  catch(std::runtime_error& e) 
  {
    std::cout << e.what() << std::endl;
    return false;
  } 

  
  forward_model = model;
  

  //any required, additional priors that are not part of the forward model
  setAdditionalPriors();


  //print the prior list to terminal
  std::cout << "\n" << "List of priors: \n";

  for (auto & i : priors)
    std::cout << i->priorName() << "\t" << i->parameterName() << "\n";

  std::cout << "\n";


  //Configure Multinest
  MultinestParameter param(config);

  size_t nb_priors = priors.size();

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


  return true;
}



//the destructor
//deletes all remaining data on the GPU
Retrieval::~Retrieval()
{
  if (config->use_gpu)
  {
    deleteFromDevice(observation_data_gpu);
    deleteFromDevice(observation_error_gpu);

    deleteFromDevice(band_sizes_gpu);
    deleteFromDevice(band_indices_gpu);
    deleteFromDevice(band_start_index_gpu);
  }

}


}
