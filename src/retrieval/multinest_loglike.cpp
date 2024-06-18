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
#include <algorithm>

#include "retrieval.h"

#include "multinest_parameter.h"
#include "priors.h"
#include "../forward_model/forward_model.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../observations/observations.h"
#include "../additional/physical_const.h"


namespace helios{


// Input arguments
// ndim 						= dimensionality (total number of free parameters) of the problem
// npars 						= total number of free plus derived parameters
// context						void pointer, any additional information
//
// Input/Output arguments
// Cube[npars] 						= on entry has the ndim parameters in unit-hypercube
//	 						on exit, the physical parameters plus copy any derived parameters you want to store with the free parameters
//
// Output arguments
// lnew 						= loglikelihood
void Retrieval::multinestLogLike(double *cube, int &ndim, int &nb_param, double &lnew, void *context)
{
  //context contains a pointer to the retrieval object
  //we now recast it accordingly to use the retrieval object here
  Retrieval *retrieval_ptr = static_cast<Retrieval*>(context);


  //convert the normalised cube values into the real parameter values
  //write back the parameter values into the cube for MultiNest output
  std::vector<double> parameter(nb_param, 0.0);

  for (size_t i=0; i<parameter.size(); ++i)
  {
    if (retrieval_ptr->priors.distributions[i]->distributionType() == "Linked prior")
      parameter[i] = retrieval_ptr->priors.distributions[i]->parameterValue(cube[retrieval_ptr->priors.prior_links[i]]);
    else
      parameter[i] = retrieval_ptr->priors.distributions[i]->parameterValue(cube[i]);
  }


  for (size_t i=0; i<parameter.size(); ++i)
    cube[i] = parameter[i];


  if (retrieval_ptr->config->multinest_print_iter_values)
  {
    std::cout << "model ";
    
    for (auto & i : parameter) std::cout << i << "   ";

    std::cout << "\n";
  }


  //run the forward model with the parameter set to obtain a high-res model spectrum
  std::vector<double> model_spectrum(retrieval_ptr->spectral_grid.nbSpectralPoints(), 0.0);
  std::vector<double> model_spectrum_bands(retrieval_ptr->nb_observation_points, 0.0);

  bool neglect = retrieval_ptr->forward_model->calcModel(parameter, model_spectrum, model_spectrum_bands);


  double error_inflation = 0;

  if (retrieval_ptr->config->use_error_inflation)
    error_inflation = std::pow(10, parameter.back());


  double log_like = 0;

  for (size_t i=0; i<retrieval_ptr->nb_observation_points; ++i)
  {
    //Eq. 22 from Paper I
    const double error_square = 
      retrieval_ptr->observation_error[i]*retrieval_ptr->observation_error[i] + error_inflation;
    
    //Eq. 23 from Paper I
    log_like += 
      (- 0.5 * std::log(error_square* 2.0 * constants::pi)
       - 0.5 * (retrieval_ptr->observation_data[i] - model_spectrum_bands[i])
             *(retrieval_ptr->observation_data[i] - model_spectrum_bands[i]) / error_square)
       * retrieval_ptr->observation_likelihood_weight[i];
  }

  
  //if the forward model tells us to neglect the current set of parameters,
  //set the likelihood to a low value
  if (neglect == true) log_like = -1e30;


  if (retrieval_ptr->config->multinest_print_iter_values)
    std::cout << log_like << "\n";


  lnew = log_like;
}





void Retrieval::multinestLogLikeGPU(double *cube, int &nb_dim, int &nb_param, double &lnew, void *context)
{
  //context contains a pointer to the retrieval object
  //we now recast it accordingly to use the retrieval object here
  Retrieval *retrieval_ptr = static_cast<Retrieval*>(context);


  if (stop_model)
  {
    std::string restart_script = "sbatch " + retrieval_ptr->config->retrieval_folder_path + "job_script_re.sh";
    std::cout << "Starting restart batch script: " << restart_script << "\n";
    std::system(restart_script.c_str());
    
    exit(0);
  }


  //convert the normalised cube values into the real parameter values
  //write back the parameter values into the cube for multinest output
  std::vector<double> parameter(nb_param, 0.0);

  for (size_t i=0; i<parameter.size(); ++i)
  {
    if (retrieval_ptr->priors.distributions[i]->distributionType() == "Linked prior")
      parameter[i] = retrieval_ptr->priors.distributions[i]->parameterValue(cube[retrieval_ptr->priors.prior_links[i]]);
    else
      parameter[i] = retrieval_ptr->priors.distributions[i]->parameterValue(cube[i]);
  }


  for (size_t i=0; i<parameter.size(); ++i)
    cube[i] = parameter[i];


  if (retrieval_ptr->config->multinest_print_iter_values)
  {
    std::cout << "model ";

    for (auto & i : parameter) std::cout << i << "   ";

    std::cout << "\n";
  }
  
  //pointer to the spectrum on the GPU
  double* model_spectrum_bands = nullptr;
  
  size_t nb_points = retrieval_ptr->spectral_grid.nbSpectralPoints();
  allocateOnDevice(model_spectrum_bands, retrieval_ptr->nb_observation_points);
  intializeOnDevice(retrieval_ptr->model_spectrum_gpu, nb_points);

  //call the forward model
  bool neglect = retrieval_ptr->forward_model->calcModelGPU(parameter, retrieval_ptr->model_spectrum_gpu, model_spectrum_bands);


  double error_inflation = 0;

  if (retrieval_ptr->config->use_error_inflation)
    error_inflation = std::pow(10, parameter.back());


  double new_log_like = retrieval_ptr->logLikeDev(model_spectrum_bands, error_inflation);

  //if the forward model tells us to neglect the current set of parameters,
  //set the likelihood to a low value
  if (neglect == true) new_log_like = -1e30;


  deleteFromDevice(model_spectrum_bands);


  if (retrieval_ptr->config->multinest_print_iter_values)
    std::cout << new_log_like << "\n";


  lnew = new_log_like;
}




//The dumper routine from MultiNest, will be called every updInt*10 iterations
//At the moment it's empty
//If you need to print/test stuff during a MultiNest run, add it here
void Retrieval::multinestDumper(
  int &nSamples,
  int &nlive,
  int &nPar,
  double **physLive,
  double **posterior,
  double **paramConstr,
  double &maxLogLike,
  double &logZ,
  double &INSlogZ,
  double &logZerr,
  void *context)
{
  
}



}
