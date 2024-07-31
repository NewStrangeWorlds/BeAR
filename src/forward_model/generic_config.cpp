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
#include <vector>
#include <sstream>
#include <algorithm>

#include "generic_config.h"
#include "../additional/exceptions.h"


namespace bear{


void GenericConfig::readAtmosphereConfig(
  std::fstream& file,
  size_t& nb_grid_points,
  std::vector<double>& pressure_boundaries)
{
  std::string line;
  std::getline(file, line); //parameter header

  std::getline(file, line);
  std::istringstream input(line);

  input >> nb_grid_points;

  std::getline(file, line); 
  std::getline(file, line); //parameter header

  std::getline(file, line);
  input.str(line); input.clear();
  
  pressure_boundaries = {0., 0.};

  input >> pressure_boundaries[0];

  std::getline(file, line); 
  std::getline(file, line); //parameter header

  std::getline(file, line);
  input.str(line); input.clear();

  input >> pressure_boundaries[1];

  std::cout << "- Atmosphere levels: " << nb_grid_points << "\n";
  std::cout << "- Pressure boundaries: ";
  for (auto & i : pressure_boundaries) std::cout << i << " ";
  std::cout << "\n";

  std::getline(file, line); //empty separator line
}


void GenericConfig::readChemistryConfig(
  std::fstream& file,
  std::vector<std::string>& models,
  std::vector<std::vector<std::string>>& parameters)
{
  std::string line;
  std::getline(file, line); //parameter header

  while (std::getline(file, line) && line.size() != 0)
  {
    std::istringstream input(line);
    
    std::string chem_model;
    input >> chem_model;

    std::cout << "- Chemistry model: " << chem_model << "\n";
    
    models.push_back(chem_model);
    parameters.resize(parameters.size()+1);

    std::string param;

    while (input >> param)
      parameters.back().push_back(param);
  }
}



void GenericConfig::readCloudConfig(
  std::fstream& file,
  std::vector<std::string>& models,
  std::vector<std::vector<std::string>>& parameters)
{
  std::string line;
  std::getline(file, line); //parameter header

  while (std::getline(file, line) && line.size() != 0)
  { 
    std::istringstream input(line);
    
    std::string model;
    input >> model;

    std::cout << "- Cloud model: " << model << "\n";

    if (model == "None" || model == "none" || model == "No" || model == "no" || model == "N" || model == "n")
      continue;
    
    models.push_back(model);
    parameters.resize(parameters.size()+1);

    std::string param;

    while (input >> param)
      parameters.back().push_back(param);
  }
}


void GenericConfig::readTemperatureConfig(
  std::fstream& file,
  std::string& model,
  std::vector<std::string>& parameters)
{
  std::string line;
  std::getline(file, line); //parameter header
  
  std::getline(file, line);
  std::istringstream input(line);
  
  std::string param = "";

  input >> model;

  while (input >> param)
    parameters.push_back(param);

  std::cout << "- Temperature profile: " << model;

  for (auto & i : parameters) std::cout << "  " << i;
  std::cout << "\n";

  std::getline(file, line); //empty separator line
}


void GenericConfig::readModuleConfig(
  std::fstream& file,
  std::vector<std::string>& modules,
  std::vector<std::vector<std::string>>& parameters)
{
  std::string line;
  std::getline(file, line); //parameter header

  while (std::getline(file, line) && line.size() != 0)
  { 
    std::istringstream input(line);
    
    std::string model;
    input >> model;

    std::cout << "- Optional modules: " << model << "\n";

    if (model == "None" || model == "none" || model == "No" || model == "no" || model == "N" || model == "n")
      continue;
    
    modules.push_back(model);
    parameters.resize(parameters.size()+1);

    std::string param;

    while (input >> param)
      parameters.back().push_back(param);
  }
}


void GenericConfig::readOpacityConfig(
  std::fstream& file,
  std::vector<std::string>& species_symbol,
  std::vector<std::string>& species_folder)
{
  std::string line;
  std::getline(file, line); //parameter header
  
  while(std::getline(file, line))
  {
    std::istringstream input(line);
    std::string species, folder;

    input >> species >> folder;
    
    if (species.length() > 0 && folder.length() > 0)
    {
      species_symbol.push_back(species);
      species_folder.push_back(folder);
    }
  }

  std::cout << "- Opacity species:\n";
  for (size_t i=0; i<species_symbol.size(); ++i)
    std::cout << "   species " << species_symbol[i] << "\t folder: " << species_folder[i] << "\n"; 
  
  std::cout << "\n";
}


bool GenericConfig::readBooleanParameter(
  std::fstream& file,
  const std::string& parameter_name)
{
  std::string line;
  std::getline(file, line);  //parameter header
  
  std::getline(file, line);
  std::istringstream input(line);
  
  std::string param = "";

  input >> param;

  std::getline(file, line); //empty separator line

  if (param == "Yes" || param == "yes" || param == "Y" || param == "y")
  {
    std::cout << "  - " << parameter_name << ": yes\n";
    return true;
  }

  if (param == "No" || param == "no" || param == "N" || param == "n")
  {
    std::cout << "  - " << parameter_name << ": no\n";
    return false;
  }
  
  std::string error_message = "Boolean parameter value for " + parameter_name + " in config file is invalid\n";
  throw InvalidInput(std::string ("GenericConfig::readBooleanParameter"), error_message);

  return false;
}


std::string GenericConfig::readParameter(
  std::fstream& file,
  const std::string& parameter_name,
  const std::vector<std::string>& allowed_values)
{
  std::string line;
  std::getline(file, line);  //parameter header
  
  std::getline(file, line);
  std::istringstream input(line);
  
  std::string param = "";

  input >> param;

  std::getline(file, line); //empty separator line

  auto it = std::find(allowed_values.begin(), allowed_values.end(), param);

  if (it == allowed_values.end())
  {
    std::string error_message = "Parameter value " + param + " for " + parameter_name + " in config file is invalid\n";
    throw InvalidInput(std::string ("GenericConfig::readParameter"), error_message);
  }
  else
    return param;
}



std::vector<chemical_species_id> GenericConfig::readChemicalSpecies(
  std::fstream& file,
  const std::string& parameter_name)
{
  std::string line;
  std::getline(file, line);  //parameter header
  
  std::getline(file, line);
  std::istringstream input(line);
  
  std::string species;
  std::vector<chemical_species_id> species_to_save;

  while (input >> species)
  {
    for (size_t j=0; j<constants::species_data.size(); ++j)
    {
      if (constants::species_data[j].symbol == species)
      {
        species_to_save.push_back(constants::species_data[j].id); 
        break;
      }
    }
  }

  std::getline(file, line); //empty separator line

  std::cout << "  - " + parameter_name + ": ";

  for (auto & i : species_to_save)
    std::cout << constants::species_data[i].symbol << " ";
  std::cout << "\n";

  return species_to_save;
}


}

