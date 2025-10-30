/*
* This file is part of the BeAR code (https://github.com/newstrangeworlds/BeAR).
* Copyright (C) 2025 Daniel Kitzmann
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



void Priors::init(
  const std::string& folder_path, 
  const size_t nb_total_param)
{
  const std::string file_name = folder_path + "priors.config";

  std::vector<std::string> prior_type; 
  std::vector<std::string> prior_description; 
  std::vector<std::vector<double>> prior_parameter;
  std::vector<std::string> prior_unit;

  readConfigFile(file_name, prior_type, prior_description, prior_parameter, prior_unit);

  if (prior_type.size() != nb_total_param)
  {
    std::string error_message = "Found " 
      + std::to_string(prior_type.size()) 
      + " priors in priors.config but expected " 
      + std::to_string(nb_total_param) + "\n";
    
    throw InvalidInput(std::string ("Priors::init"), error_message);
  }
  
  std::vector<PriorConfig> priors_config;
  
  for (size_t i=0; i<prior_type.size(); ++i)
    priors_config.push_back(PriorConfig(prior_type[i], prior_description[i], prior_parameter[i], prior_unit[i]));

  add(priors_config);
}



void Priors::init(
  const std::vector<PriorConfig>& priors_config,
  const size_t nb_total_param)
{
  if (priors_config.size() != nb_total_param)
  {
    std::string error_message = "Found " 
      + std::to_string(priors_config.size()) 
      + " priors config but expected " 
      + std::to_string(nb_total_param) + "\n";
    
    throw InvalidInput(std::string ("Priors::init"), error_message);
  }

  add(priors_config);
}



void Priors::printInfo()
{
  std::cout << "\n" << "List of priors: \n";

  for (auto & i : distributions)
    i->printInfo();

  std::cout << "\n";
}



void Priors::add(
  const std::vector<PriorConfig>& priors_config)
{
  std::vector<std::string> type(priors_config.size(), "");
  std::vector<std::string> description(priors_config.size(), "");
  std::vector<std::vector<double>> parameter(priors_config.size(), std::vector<double>());
  std::vector<std::string> unit(priors_config.size(), "");
  
  for (size_t i=0; i<priors_config.size(); ++i)
  {
    type[i] = priors_config[i].type;
    description[i] = priors_config[i].description;

    for (size_t j=0; j<priors_config[i].parameter.size(); ++j)
      parameter[i].push_back(priors_config[i].parameter[j]);
    
    unit[i] = priors_config[i].unit;
  }
  
  for (size_t i=0; i<priors_config.size(); ++i)
    addSingle(type[i], description[i], parameter[i], unit[i]);

  setupLinkedPriors(type, description, parameter);
}



void Priors::addSingle(
  const std::string& type, 
  const std::string& description, 
  const std::vector<double>& parameter,
  const std::string& unit)
{
  if (type == "uniform")
  {
    if (parameter.size() != 2)
    {
      std::string error_message = "uniform prior " + description + " requires two parameters with an optional unit!\n";
      throw InvalidInput(std::string ("pirors.config"), error_message);
    }

    if (unit == "")
    {
      UniformPrior* uniform_prior = new UniformPrior(
        description, parameter[0], parameter[1], "");
      distributions.push_back(uniform_prior);
    }
    else
    {
      UniformPrior* uniform_prior = new UniformPrior(
        description, parameter[0], parameter[1], unit);
      distributions.push_back(uniform_prior);
    }

    return;
  }


  if (type == "log_uniform")
  {
    if (parameter.size() != 2)
    {
      std::string error_message = "log uniform prior " + description + " requires two parameters with an optional unit!\n";
      throw InvalidInput(std::string ("pirors.config"), error_message);
    }

    if (unit == "")
    {
      LogUniformPrior* uniform_prior = new LogUniformPrior(
        description, parameter[0], parameter[1], "");
      distributions.push_back(uniform_prior);
    }
    else
    {
      LogUniformPrior* uniform_prior = new LogUniformPrior(
        description, parameter[0], parameter[1], unit);
      distributions.push_back(uniform_prior);
    }

    return;
  }


  if (type == "gaussian")
  {
    if (parameter.size() != 2)
    {
      std::string error_message = "Gaussian prior " + description + " requires two parameters with an optional unit!\n";
      throw InvalidInput(std::string ("pirors.config"), error_message);
    }

    if (unit == "")
    {
      GaussianPrior* gaussian_prior = new GaussianPrior(
        description, parameter[0], parameter[1], "");
      distributions.push_back(gaussian_prior);
    }
    else
    {
      GaussianPrior* gaussian_prior = new GaussianPrior(
        description, parameter[0], parameter[1], unit);
      distributions.push_back(gaussian_prior);
    }

    return;
  }


  if (type == "delta")
  {
    if (parameter.size() != 1)
    {
      std::string error_message = "Delta distribution prior " + description + " requires one parameter with an optional unit!\n";
      throw InvalidInput(std::string ("pirors.config"), error_message);
    }

    if (unit == "")
    {
      DeltaPrior* delta_prior = new DeltaPrior(
        description, parameter[0], "");
      distributions.push_back(delta_prior);
    }
    else
    {
      DeltaPrior* delta_prior = new DeltaPrior(
        description, parameter[0], unit);
      distributions.push_back(delta_prior);
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
    
    DeltaPrior* delta_prior = new DeltaPrior(
      description, parameter[0], "");
    distributions.push_back(delta_prior);

    return;
  }

  std::string error_message = "Prior type " + type + " for " + description + " unknown!\n";
  throw InvalidInput(std::string ("Retrieval::setPrior"), error_message);
}




void Priors::setupLinkedPriors(
  const std::vector<std::string>& type,
  const std::vector<std::string>& description,
  const std::vector<std::vector<double>>& parameter)
{
  prior_links.assign(type.size(), 0);

  for (size_t i=0; i<type.size(); ++i)
  {
    if (type[i] == "linked")
    {
      //delete the place holder
      delete(distributions[i]);
      
      size_t prior_index = static_cast<int>(parameter[i][0]) - 1;

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



void Priors::readConfigFile(
  const std::string& file_path, 
  std::vector<std::string>& prior_type, 
  std::vector<std::string>& prior_description, 
  std::vector<std::vector<double>>& prior_parameter,
  std::vector<std::string>& prior_unit)
{

  std::fstream file;
  file.open(file_path.c_str(), std::ios::in);

  if (file.fail())  
    throw FileNotFound(std::string ("ForwardModel::readPriorConfigFile"), file_path);


  auto is_number = [](const std::string& s){
    std::istringstream iss(s);
    double d;
    return iss >> std::noskipws >> d && iss.eof();};


  std::string line;

  while (std::getline(file, line))
  {
    std::istringstream input(line);

    std::string type, description;
    std::vector<std::string> parameter;
    std::string unit = "";

    input >> type >> description;

    std::string single_parameter;

    while (input >> single_parameter)
      parameter.push_back(single_parameter);

    std::vector<double> double_parameter;

    if (is_number(parameter.back()) == true)
    {
      for (size_t i=0; i<parameter.size(); ++i)
      {
        if (is_number(parameter[i]) == false)
        {
          std::string error_message = "Prior value " + parameter[i] + " cannot be converted into a number!\n";
          throw InvalidInput(std::string ("pirors.config"), error_message);
        }
        
        double_parameter.push_back(std::stod(parameter[i]));
      }
        
    }
    else
    {
      for (size_t i=0; i<parameter.size()-1; ++i)
      {
        if (is_number(parameter[i]) == false)
        {
          std::string error_message = "Prior value " + parameter[i] + " cannot be converted into a number!\n";
          throw InvalidInput(std::string ("pirors.config"), error_message);
        }

        double_parameter.push_back(std::stod(parameter[i]));
      }
      
      unit = parameter.back();
    }

    prior_type.push_back(type);
    prior_description.push_back(description);
    prior_parameter.push_back(double_parameter);
    prior_unit.push_back(unit);
  }

  file.close();
}


}
