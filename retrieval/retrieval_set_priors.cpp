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



namespace helios{



//Defines the prior distributions for all parameters
void Retrieval::setPriors()
{
  UniformPrior* uniform_prior = nullptr;
  LogUniformPrior* log_uniform_prior = nullptr;
  GaussianPrior* gaussian_prior = nullptr;


  uniform_prior = new UniformPrior("log g", 3.5, 6.0);
  priors.push_back(uniform_prior);

  uniform_prior = new UniformPrior("scaling factor", 0.1, 5.0);
  priors.push_back(uniform_prior);
  
  gaussian_prior = new GaussianPrior("distance", 3.6389, 0.0033);
  priors.push_back(gaussian_prior);


  log_uniform_prior = new LogUniformPrior("MR H2O", 1e-12, 0.1);
  priors.push_back(log_uniform_prior);

  log_uniform_prior = new LogUniformPrior("MR CH4", 1e-12, 0.1);
  priors.push_back(log_uniform_prior);

  log_uniform_prior = new LogUniformPrior("MR NH3", 1e-12, 0.1);
  priors.push_back(log_uniform_prior);

  log_uniform_prior = new LogUniformPrior("MR CO2", 1e-12, 0.1);
  priors.push_back(log_uniform_prior);

  log_uniform_prior = new LogUniformPrior("MR CO", 1e-12, 0.1);
  priors.push_back(log_uniform_prior);

  log_uniform_prior = new LogUniformPrior("MR H2S", 1e-12, 0.1);
  priors.push_back(log_uniform_prior);

  log_uniform_prior = new LogUniformPrior("MR Na", 1e-12, 0.001);
  priors.push_back(log_uniform_prior);



  uniform_prior = new UniformPrior("temperature 1", 3000, 1000);
  priors.push_back(uniform_prior);

  uniform_prior = new UniformPrior("temperature 2", 0.3, 0.95);
  priors.push_back(uniform_prior);

  uniform_prior = new UniformPrior("temperature 2", 0.3, 0.95);
  priors.push_back(uniform_prior);

  uniform_prior = new UniformPrior("temperature 2", 0.4, 0.95);
  priors.push_back(uniform_prior);

  uniform_prior = new UniformPrior("temperature 3", 0.5, 0.95);
  priors.push_back(uniform_prior);

  uniform_prior = new UniformPrior("temperature 4", 0.5, 0.95);
  priors.push_back(uniform_prior);

  uniform_prior = new UniformPrior("temperature 5", 0.5, 0.95);
  priors.push_back(uniform_prior);


  /*log_uniform_prior = new LogUniformPrior("cloud top", 1e-2, 50);
  priors.push_back(log_uniform_prior);

  log_uniform_prior = new LogUniformPrior("cloud bottom", 1, 10);
  priors.push_back(log_uniform_prior);

  log_uniform_prior = new LogUniformPrior("cloud optical depth", 1e-5, 20);
  priors.push_back(log_uniform_prior);*/



  //this creates the prior distribution for the error exponent
  //first, we need to find the minimum and maximum values of the observational data errors
  std::vector<double>::iterator it = std::min_element(std::begin(observation_error), std::end(observation_error));
  const double error_min = std::log10(0.1 * *it * *it);

  it = std::max_element(std::begin(observation_error), std::end(observation_error));
  const double error_max = std::log10(100.0 * *it * *it);

  //then we create the corresponding prior
  uniform_prior = new UniformPrior("error exponent", error_min, error_max);
  priors.push_back(uniform_prior);


  //print the prior list to terminal
  std::cout << "\n" << "List of priors: \n";

  for (size_t i=0; i<priors.size(); ++i)
    std::cout << "prior: " << priors[i]->priorName() << "\t" << priors[i]->parameterName() << "\n";

  std::cout << "\n";
}




}
