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


namespace bear{


//convert the normalised cube values into the real parameter values
//write back the parameter values into the cube for MultiNest output
std::pair<std::vector<double>, std::vector<double>> Retrieval::convertCubeParameters(
  std::vector<double>& cube)
{
  std::vector<double> parameter;
  std::vector<double> physical_parameter;

  convertHypercubeParameters(
    cube.data(),
    cube.size(),
    parameter,
    physical_parameter);
  
  return std::make_pair(parameter, physical_parameter);
}


//convert the normalised cube values into the real parameter values
//write back the parameter values into the cube for MultiNest output
void Retrieval::convertHypercubeParameters(
  double *cube,
  const size_t nb_param,
  std::vector<double>& parameter,
  std::vector<double>& physical_parameter)
{
  parameter.assign(nb_param, 0.0);
  physical_parameter.assign(nb_param, 0.0);

  for (size_t i=0; i<parameter.size(); ++i)
  {
    if (priors.distributions[i]->distributionType() == "Linked prior")
    {
      parameter[i] = priors.distributions[i]->parameterValue(
        cube[priors.prior_links[i]]);
      physical_parameter[i] = priors.distributions[i]->parameterPhysicalValue(
        cube[priors.prior_links[i]]);
    }
    else
    {
      parameter[i] = priors.distributions[i]->parameterValue(cube[i]);
      physical_parameter[i] = priors.distributions[i]->parameterPhysicalValue(cube[i]);
    }
  }

  for (size_t i=0; i<parameter.size(); ++i)
    cube[i] = parameter[i];

  if (config->multinest_print_iter_values)
  {
    std::cout << "model ";
    
    for (auto & i : parameter) std::cout << i << "   ";
    std::cout << "\n";

    for (auto & i : physical_parameter) std::cout << i << "   ";
    std::cout << "\n";
  }
}


std::vector<double> Retrieval::convertToPhysicalParameters(
  const std::vector<double>& parameters)
{
  std::vector<double> physical_parameters(parameters.size(), 0.0);
  
  if (parameters.size() != priors.number())
  {
    std::string error_message = 
      "Number of posterior parameters not equal to the number of free parameters of the forward model.\n";
    throw InvalidInput(std::string ("Retrieval::convertToPhysicalParameters"), error_message);
  }

  for (size_t i=0; i<parameters.size(); ++i)
    physical_parameters[i] = priors.distributions[i]->applyParameterUnit(parameters[i]);

  return physical_parameters;
}



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
void Retrieval::multinestLogLike(
  double *cube, int &ndim, int &nb_param, double &lnew, void *context)
{
  //context contains a pointer to the retrieval object
  //we now recast it accordingly to use the retrieval object here
  Retrieval *retrieval_ptr = static_cast<Retrieval*>(context);
  
  std::vector<double> parameters {};
  std::vector<double> physical_parameters {};

  retrieval_ptr->convertHypercubeParameters(
    cube, 
    nb_param, 
    parameters, 
    physical_parameters);

  lnew = retrieval_ptr->logLikelihood(physical_parameters);
}



void Retrieval::multinestLogLikeGPU(
  double *cube, int &nb_dim, int &nb_param, double &lnew, void *context)
{
  //context contains a pointer to the retrieval object
  //we now recast it accordingly to use the retrieval object here
  Retrieval *retrieval_ptr = static_cast<Retrieval*>(context);


  if (stop_model)
  {
    std::string restart_script = "sbatch " 
      + retrieval_ptr->config->retrieval_folder_path 
      + "job_script_re.sh";
    std::cout << "Starting restart batch script: " << restart_script << "\n";
    std::system(restart_script.c_str());
    
    exit(0);
  }


  std::vector<double> parameter {};
  std::vector<double> physical_parameter {};

  retrieval_ptr->convertHypercubeParameters(
    cube, 
    nb_param, 
    parameter, 
    physical_parameter);
  
  lnew = retrieval_ptr->logLikelihoodGPU(physical_parameter);
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
