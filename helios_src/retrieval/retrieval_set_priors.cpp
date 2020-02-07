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


#include "multinest_parameter.h"
#include "prior.h"
#include "../observations/observations.h"
#include "../additional/exceptions.h"



namespace helios{



void Retrieval::setPriors(const std::vector<std::string>& type, const std::vector<std::string>& description, const std::vector<std::vector<double>>& parameter)
{
  if (description.size() != type.size() || parameter.size() != type.size())
  {
    std::string error_message = "Incorrect prior data sizes! Type, description, and parameter should have the same length.\n";
    throw ExceptionInvalidInput(std::string ("Retrieval::setPriors"), error_message);
  }


  for (size_t i=0; i<type.size(); ++i)
    setPrior(type[i], description[i], parameter[i]);
}




void Retrieval::setPrior(const std::string& type, const std::string& description, const std::vector<double>& parameter)
{
  //find the corresponding prior to the supplied "type" string
  auto it = std::find(helios::priors::prior_type_strings.begin(), helios::priors::prior_type_strings.end(), type);

  if (it == helios::priors::prior_type_strings.end())
  {
    std::string error_message = "Prior type " + type + " for " + description + " unknown!\n";
    throw ExceptionInvalidInput(std::string ("Retrieval::setPrior"), error_message);
  }


  PriorType prior_type = helios::priors::prior_types[std::distance(helios::priors::prior_type_strings.begin(), it)];


  UniformPrior* uniform_prior = nullptr;
  LogUniformPrior* log_uniform_prior = nullptr;
  GaussianPrior* gaussian_prior = nullptr;

  
  switch (prior_type)
  {
    case PriorType::uniform :     if (parameter.size() == 2) 
                                  {
                                    uniform_prior = new UniformPrior(description, parameter[0], parameter[1]);
                                    priors.push_back(uniform_prior);
                                  }
                                  else
                                  { 
                                    std::string error_message = "uniform prior " + description + " requires two parameters!\n"; 
                                    throw ExceptionInvalidInput(std::string ("Retrieval::setPrior"), error_message);
                                  }
                                  break;
    case PriorType::log_uniform : if (parameter.size() == 2)
                                  {
                                    log_uniform_prior = new LogUniformPrior(description, parameter[0], parameter[1]);
                                    priors.push_back(log_uniform_prior);
                                  }    
                                  else
                                  { 
                                    std::string error_message = "log-uniform prior " + description + " requires two parameters!\n"; 
                                    throw ExceptionInvalidInput(std::string ("Retrieval::setPrior"), error_message);
                                  }
                                  break;
    case PriorType::gaussian :    if (parameter.size() == 2)
                                  {
                                    gaussian_prior = new GaussianPrior(description, parameter[0], parameter[1]);
                                    priors.push_back(gaussian_prior);
                                  } 
                                  else
                                  { 
                                    std::string error_message = "gaussian prior " + description + " requires two parameters!\n"; 
                                    throw ExceptionInvalidInput(std::string ("Retrieval::setPrior"), error_message);
                                  }  
                                  break;
    default : std::cout << "Unknown prior type " << type << " for " << description << "\n"; //this shouldn't really happen here
  }

}



//Add additional priors that are not set by the forward model
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

    //then we create the corresponding prior
    UniformPrior* uniform_prior = new UniformPrior("error exponent", error_min, error_max);
    priors.push_back(uniform_prior);
  }

}




}
