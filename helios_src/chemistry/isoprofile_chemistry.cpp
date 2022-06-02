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


#include "isoprofile_chemistry.h"
#include "chem_species.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"

#include <algorithm>
#include <vector>



namespace helios {


IsoprofileChemistry::IsoprofileChemistry(const std::vector<std::string>& chemical_species)
{ 
  std::cout << "- Chemistry model: " << "isoprofiles" << "\n";
  std::cout << "  - Species for this model: ";

  for (auto & i : chemical_species)
    std::cout << i << "  ";
  
  std::cout << "\n";

  for (auto & i : chemical_species)
  {
    //H2 and He are used as a buffer gas and are not free variables
    if (i == "H2" || i == "He")
    {
      std::string error_message = "Chemical species " + i + " is used as a buffer gas and is not a free variable in IsoprofileChemistry. It has to be removed from the config file and as a prior!\n";
      throw ExceptionInvalidInput(std::string ("IsoprofileChemistry::IsoprofileChemistry"), error_message);
    } 

    bool species_found = false;

    for (size_t j=0; j<constants::species_data.size(); ++j)
    {
      if (constants::species_data[j].symbol == i)
      {
        species.push_back(constants::species_data[j].id); 
        species_found = true;
        break;
      } 
    }


    if (!species_found)
    {
      std::string error_message = "Chemical species " + i + " not found in the list of species in chem_species.h \n";
      throw ExceptionInvalidInput(std::string ("IsoprofileChemistry::IsoprofileChemistry"), error_message);
    }
      
  }

  /*std::cout << "Isoprofile chemistry initialised with the following species: \n";
  for (auto & i : species)
    std::cout << constants::species_data[i].symbol << "  ";
  std::cout << "\n";*/


  //check if Na is a free variable
  auto it = std::find(species.begin(), species.end(), _Na);
  
  if (it == species.end()) 
    sodium_free_parameter = false;
  else
    sodium_free_parameter = true;

  
  nb_parameters = chemical_species.size();
}




bool IsoprofileChemistry::calcChemicalComposition(const std::vector<double>& parameters, const std::vector<double>& temperature, const std::vector<double>& pressure,
                                                  std::vector<std::vector<double>>& number_densities, std::vector<double>& mean_molecular_weight)
{
  const double solar_h2_he = solar_h2 + solar_he;
  const double solar_na_k = 16.2181;
  
  const double solar_h_he = solar_he + 1.0;
  const double epsilon_h = 1.0 / (solar_h_he);

  bool neglect_model = false;

  for (size_t i=0; i<number_densities.size(); ++i)
  { 
    number_densities[i][_TOTAL] = pressure[i] * 1.e6 / constants::boltzmann_k / temperature[i];

    for (size_t j=0; j<species.size(); ++j)
      number_densities[i][species[j]] = number_densities[i][_TOTAL] * parameters[j];

    //for (size_t j=0; j<species.size(); ++j)
      //number_densities[i][species[j]] *= parameters[j];

    if (!sodium_free_parameter)
      number_densities[i][_Na] = number_densities[i][_K] * solar_na_k;


    //the mixing ratio of the two background gases H2 and He
    double mixing_ratio_background = 1.0;

    for (auto & j : constants::species_data)
      if (j.id != _H2 && j.id != _He && j.id != _TOTAL && j.id != _H) 
        mixing_ratio_background -= number_densities[i][j.id]/number_densities[i][_TOTAL];


    //if we have a negative mixing ratio of the background gas
    //(i.e. the sum of all other species is larger than 1)
    //neglect this model
    if (mixing_ratio_background < 0)
    {
      neglect_model = true;
      mixing_ratio_background = 0;
    }
    
    //const double mixing_ratio_h = number_densities[i][_H]/number_densities[i][_TOTAL];
    //const double mixing_ratio_h2 = epsilon_h / (2.0 - epsilon_h) * (mixing_ratio_background - mixing_ratio_h);
    //const double mixing_ratio_he = mixing_ratio_background - mixing_ratio_h - mixing_ratio_h2;

    const double mixing_ratio_h2 = mixing_ratio_background * solar_h2 / solar_h2_he;
    const double mixing_ratio_he = mixing_ratio_background * solar_he / solar_h2_he;

    number_densities[i][_H2] = number_densities[i][_TOTAL] * mixing_ratio_h2;
    number_densities[i][_He] = number_densities[i][_TOTAL] * mixing_ratio_he;
  }


  //calculate the mean molecular weight
  //note that we sum over *all* species, not just the ones that were included in this chemistry
  for (size_t i=0; i<number_densities.size(); ++i)
  {
    double mu = 0;

    for (auto & j : constants::species_data)
      mu += number_densities[i][j.id]/number_densities[i][_TOTAL] * j.molecular_weight;

    mean_molecular_weight[i] = mu;
  }


  return neglect_model;
}



}

