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


#include "secondary_eclipse.h"


#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>


#include "../../retrieval/priors.h"
#include "../../additional/exceptions.h"


namespace helios{


//set the model priors
void SecondaryEclipseModel::setPriors(Priors* priors)
{
  const std::string file_name = config->retrieval_folder_path + "priors.config";

  std::vector<std::string> prior_type; 
  std::vector<std::string> prior_description; 
  std::vector<std::vector<double>> prior_parameter;


  readPriorConfigFile(file_name, prior_type, prior_description, prior_parameter);


  //check if we have the correct number of piors
  if (prior_type.size() != nb_total_param())
  {
    std::string error_message = "Found " + std::to_string(prior_type.size()) + " priors in priors.config but expected " + std::to_string(nb_total_param()) + "\n";
    throw ExceptionInvalidInput(std::string ("SecondaryEclipseModel::setPriors"), error_message);
  }


  priors->add(prior_type, prior_description, prior_parameter);
}




void SecondaryEclipseModel::readPriorConfigFile(const std::string& file_path, std::vector<std::string>& prior_type, 
                                                                        std::vector<std::string>& prior_description, 
                                                                        std::vector<std::vector<double>>& prior_parameter)
{
  std::fstream file;
  file.open(file_path.c_str(), std::ios::in);


  if (file.fail())  
    throw ExceptionFileNotFound(std::string ("SecondaryEclipseModel::readPriorConfigFile"), file_path);


  std::string line;

  while (std::getline(file, line))
  {
    std::istringstream input(line);

    std::string type, description;
    std::vector<double> parameter;

    input >> type >> description;

    double single_parameter;

    while (input >> single_parameter)
      parameter.push_back(single_parameter);


    prior_type.push_back(type);
    prior_description.push_back(description);
    prior_parameter.push_back(parameter);
  }


  file.close();
}



}

