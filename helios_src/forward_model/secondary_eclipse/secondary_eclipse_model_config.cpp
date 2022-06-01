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


#include "../../retrieval/retrieval.h"
#include "../../retrieval/prior.h"
#include "../../additional/exceptions.h"



namespace helios{


SecondaryEclipseConfig::SecondaryEclipseConfig (const std::string& folder_path)
{
  const std::string config_file_name = folder_path + "forward_model.config";

  readConfigFile(config_file_name);
}



void SecondaryEclipseConfig::readConfigFile(const std::string& file_name)
{
  std::fstream file;
  file.open(file_name.c_str(), std::ios::in);

  
  if (file.fail())  
    throw ExceptionFileNotFound(std::string ("SecondaryEclipseConfig::readConfigFile"), file_name);

  
  std::string line;
  std::string input;


  std::getline(file, line);
  
  file >> nb_grid_points >> line;
  std::cout << "- Atmosphere levels: " << nb_grid_points << "\n";


  std::getline(file, line);

  file >> atmos_top_pressure >> line;
  std::cout << "- Top of atmosphere pressure: " << atmos_top_pressure << "\n";


  std::getline(file, line);

  file >> atmos_bottom_pressure >> line;
  std::cout << "- Bottom of atmosphere pressure: " << atmos_bottom_pressure << "\n";


  std::getline(file, line);

  file >> nb_temperature_elements >> line;
  std::cout << "- Number of Temperature elements: " << nb_temperature_elements << "\n";
  

  atmos_boundaries[0] = atmos_top_pressure;
  atmos_boundaries[1] = atmos_bottom_pressure;


  std::getline(file, line);

  file >> temperature_poly_degree >> line;
  std::cout << "- Temperature polynomial degree: " << temperature_poly_degree << "\n";
  
  
  std::getline(file, line);

  file >> stellar_spectrum_file >> line;
  std::cout << "- Stellar spectrum file: " << stellar_spectrum_file << "\n";


  std::getline(file, line);

  file >> input >> line;
  if (input == "Y" || input == "Yes" || input == "1") use_cloud_layer = true;
  std::cout << "- Use Cloud Layer: " << use_cloud_layer << "\n";


  std::getline(file, line);

  file >> input >> line;
  
  if (input == "scm") radiative_transfer_model = 0;
  if (input == "disort") radiative_transfer_model = 1;
  if (input != "disort" && input != "scm")
  {
    std::string error_message = "Radiative transfer model " + input + " in forward_model.config unknown!\n";
    throw ExceptionInvalidInput(std::string ("SecondaryEclipseConfig::readConfigFile"), error_message);
  }

  std::cout << "- Radiative transfer model: " << radiative_transfer_model << "\n";


  readChemistryConfig(file);

  readOpacityConfig(file);
  

  file.close();
}



void SecondaryEclipseConfig::readChemistryConfig(std::fstream& file)
{
  std::string line;
  std::getline(file, line);  
  

  while (std::getline(file, line) && line.size() != 0)
  { 
    std::istringstream input(line);

    std::string chem_model;
    input >> chem_model;

    if (chem_model == "iso")
    {
      chemistry_model.push_back(0);

      chemistry_parameters.resize(chemistry_parameters.size()+1);
    
      std::string species;

      while (input >> species)
        chemistry_parameters.back().push_back(species);

      std::cout << "- Chemistry model: " << "isoprofiles" << "\n";
      std::cout << "  - Species for this model: ";
    
      for (auto & i : chemistry_parameters.back())
        std::cout << i << "  ";

      std::cout << "\n";
    }


    if (chem_model == "eq")
    {
      chemistry_model.push_back(1);

      chemistry_parameters.push_back(std::vector<std::string>(1, ""));

      input >> chemistry_parameters.back()[0];
    
      std::cout << "- Chemistry model: " << "FastChem" << "\n";
      std::cout << "  - Parameter file: " << chemistry_parameters.back()[0] << "\n";
    }


    if (chem_model == "free")
    {
      chemistry_model.push_back(2);

      chemistry_parameters.push_back(std::vector<std::string>(3, ""));

      input >> chemistry_parameters.back()[0] >> chemistry_parameters.back()[1] >> chemistry_parameters.back()[2];
    
      std::cout << "- Chemistry model: " << "free chemistry" << "\n";
      std::cout << "  - Species for this model: " << chemistry_parameters.back()[0] << "\n";
      std::cout << "  - number of elements: " << chemistry_parameters.back()[1] << "\n";
      std::cout << "  - polynomial degree: " << chemistry_parameters.back()[2] << "\n";
    }


    if (chem_model != "eq" && chem_model != "iso" && chem_model != "free")
    {
      std::string error_message = "Chemistry model " + chem_model + " in forward_model.config unknown!\n";
      throw ExceptionInvalidInput(std::string ("SecondaryEclipseConfig::readConfigFile"), error_message);
    }

  }
 
}



void SecondaryEclipseConfig::readOpacityConfig(std::fstream& file)
{
  std::string line;
  std::getline(file, line);
  
  
  while(std::getline(file, line))
  {
    std::istringstream input(line);

    std::string species, folder;

    input >> species >> folder;
    
    if (species.length() > 0 && folder.length() > 0)
    {
      opacity_species_symbol.push_back(species);
      opacity_species_folder.push_back(folder);
    }
    
  }

  
  std::cout << "- Opacity species:\n";
  for (size_t i=0; i<opacity_species_symbol.size(); ++i)
    std::cout << "   species " << opacity_species_symbol[i] << "\t folder: " << opacity_species_folder[i] << "\n"; 
  
  
  std::cout << "\n";
}



}

