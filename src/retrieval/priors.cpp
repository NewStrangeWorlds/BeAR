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
#include <cmath>
#include <vector>
#include <algorithm>

#include "priors.h"

#include "prior_types.h"
#include "../additional/exceptions.h"


namespace bear{


Priors::~Priors()
{
   for (size_t i=0; i<distributions.size(); ++i)
    delete distributions[i];
}



void Priors::add(
  const std::vector<std::string>& type, 
  const std::vector<std::string>& description, 
  const std::vector<std::vector<std::string>>& parameter)
{
  if (description.size() != type.size() || parameter.size() != type.size())
  {
    std::string error_message = "Incorrect prior data sizes! Type, description, and parameter should have the same length.\n";
    throw InvalidInput(std::string ("Retrieval::setPriors"), error_message);
  }


  for (size_t i=0; i<type.size(); ++i)
    addSingle(type[i], description[i], parameter[i]);

  setupLinkedPriors(type, description, parameter);
}


void Priors::addSingle(
  const std::string& type, 
  const std::string& description, 
  const std::vector<std::string>& parameter)
{
  if (type == "uniform")
  {
    if (parameter.size() != 2 && parameter.size() != 3)
    {
      std::string error_message = "uniform prior " + description + " requires two parameters with an optional unit!\n";
      throw InvalidInput(std::string ("pirors.config"), error_message);
    }

    if (parameter.size() == 2)
    {
      UniformPrior* uniform_prior = new UniformPrior(
        description, std::stod(parameter[0]), std::stod(parameter[1]), "");
      distributions.push_back(uniform_prior);
    }
    else
    {
      UniformPrior* uniform_prior = new UniformPrior(
        description, std::stod(parameter[0]), std::stod(parameter[1]), parameter[2]);
      distributions.push_back(uniform_prior);
    }

    return;
  }


  if (type == "log_uniform")
  {
    if (parameter.size() != 2 && parameter.size() != 3)
    {
      std::string error_message = "log uniform prior " + description + " requires two parameters with an optional unit!\n";
      throw InvalidInput(std::string ("pirors.config"), error_message);
    }

    if (parameter.size() == 2)
    {
      LogUniformPrior* uniform_prior = new LogUniformPrior(
        description, std::stod(parameter[0]), std::stod(parameter[1]), "");
      distributions.push_back(uniform_prior);
    }
    else
    {
      LogUniformPrior* uniform_prior = new LogUniformPrior(
        description, std::stod(parameter[0]), std::stod(parameter[1]), parameter[2]);
      distributions.push_back(uniform_prior);
    }

    return;
  }


  if (type == "gaussian")
  {
    if (parameter.size() != 2 && parameter.size() != 3)
    {
      std::string error_message = "Gaussian prior " + description + " requires two parameters with an optional unit!\n";
      throw InvalidInput(std::string ("pirors.config"), error_message);
    }

    if (parameter.size() == 2)
    {
      GaussianPrior* gaussian_prior = new GaussianPrior(
        description, std::stod(parameter[0]), std::stod(parameter[1]), "");
      distributions.push_back(gaussian_prior);
    }
    else
    {
      GaussianPrior* uniform_prior = new GaussianPrior(
        description, std::stod(parameter[0]), std::stod(parameter[1]), parameter[2]);
      distributions.push_back(uniform_prior);
    }

    return;
  }


  if (type == "delta")
  {
    if (parameter.size() != 1 && parameter.size() != 2)
    {
      std::string error_message = "Delta distribution prior " + description + " requires one parameter with an optional unit!\n";
      throw InvalidInput(std::string ("pirors.config"), error_message);
    }

    if (parameter.size() == 1)
    {
      DeltaPrior* gaussian_prior = new DeltaPrior(
        description, std::stod(parameter[0]), "");
      distributions.push_back(gaussian_prior);
    }
    else
    {
      DeltaPrior* uniform_prior = new DeltaPrior(
        description, std::stod(parameter[0]), parameter[1]);
      distributions.push_back(uniform_prior);
    }

    return;
  }


  if (type == "linked")
  {
    if (parameter.size() != 1)
    {
      std::string error_message = "linked prior " + description + " requires one parameter!\n";
      throw InvalidInput(std::string ("pirors.config"), error_message);
    }
    
    DeltaPrior* gaussian_prior = new DeltaPrior(
      description, std::stod(parameter[0]), "");
    distributions.push_back(gaussian_prior);

    return;
  }

  std::string error_message = "Prior type " + type + " for " + description + " unknown!\n";
  throw InvalidInput(std::string ("Retrieval::setPrior"), error_message);
}




void Priors::printInfo()
{
  std::cout << "\n" << "List of priors: \n";

  for (auto & i : distributions)
    i->printInfo();

  std::cout << "\n";
}



void Priors::setupLinkedPriors(
  const std::vector<std::string>& type,
  const std::vector<std::string>& description,
  const std::vector<std::vector<std::string>>& parameter)
{
  prior_links.assign(type.size(), 0);

  for (size_t i=0; i<type.size(); ++i)
  {
    if (type[i] == "linked")
    {
      //delete the place holder
      delete(distributions[i]);
      
      size_t prior_index = std::stoi(parameter[i][0]) - 1;

      if (prior_index > type.size())
      {
        std::string error_message = "Linked prior " + description[i] + " wrong index for link!\n";
        throw InvalidInput(std::string ("pirors.config"), error_message);
      }
      
      if (type[prior_index] == "linked")
      {
        std::string error_message = "Linked prior " + description[i] + " Can not link prior to another linked prior!\n";
        throw InvalidInput(std::string ("pirors.config"), error_message);
      }
 
      LinkedPrior* linked_prior = new LinkedPrior(description[i], distributions[prior_index]);
      distributions[i] = linked_prior;
      prior_links[i] = prior_index;
    }
  }
}


}
