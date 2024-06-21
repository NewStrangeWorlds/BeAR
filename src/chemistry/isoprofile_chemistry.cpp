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


#include <algorithm>
#include <vector>

#include "isoprofile_chemistry.h"

#include "chem_species.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"


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
      throw InvalidInput(std::string ("IsoprofileChemistry::IsoprofileChemistry"), error_message);
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



bool IsoprofileChemistry::calcChemicalComposition(
  const std::vector<double>& parameters,
  const std::vector<double>& temperature,
  const std::vector<double>& pressure,
  std::vector<std::vector<double>>& number_densities,
  std::vector<double>& mean_molecular_weight)
{
  const double solar_na_k = 16.2181;

  bool neglect_model = false;

  for (size_t i=0; i<number_densities.size(); ++i)
  { 
    number_densities[i][_TOTAL] = pressure[i] * 1.e6 / constants::boltzmann_k / temperature[i];

    for (size_t j=0; j<species.size(); ++j)
      number_densities[i][species[j]] = number_densities[i][_TOTAL] * parameters[j];

    if (!sodium_free_parameter)
      number_densities[i][_Na] = number_densities[i][_K] * solar_na_k;
  }

  //calculate the mean molecular weight
  //note that we sum over *all* species, not just the ones that were included in this chemistry
  meanMolecularWeight(number_densities, mean_molecular_weight);
  
  bool mixing_ratios_ok = checkMixingRatios(number_densities);

  if (!mixing_ratios_ok)
    neglect_model = true;

  return neglect_model;
}



}

