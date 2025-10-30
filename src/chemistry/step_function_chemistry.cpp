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
#include <cmath>

#include "step_function_chemistry.h"

#include "chem_species.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"


namespace bear {


StepFunctionChemistry::StepFunctionChemistry(
  const std::string& chemical_species)
{ 
  std::cout << "- Chemistry model: " << "step-function chemistry" << "\n";
  std::cout << "  - Species for this model: " << chemical_species << "\n";
  std::cout << "\n";

  bool species_found = false;

  for (size_t j=0; j<constants::species_data.size(); ++j)
  {
    if (constants::species_data[j].symbol == chemical_species)
    {
      species.push_back(constants::species_data[j].id); 
      species_found = true;
      break;
    } 
  }

  if (!species_found)
  {
    std::string error_message = "Chemical species " + chemical_species + " not found in the list of species in chem_species.h \n";
    throw InvalidInput(std::string ("StepFunctionChemistry::StepFunctionChemistry"), error_message);
  }

  nb_parameters = 3;
}



//the parameters should contain the mixing ratios at the dof
bool StepFunctionChemistry::calcChemicalComposition(
  const std::vector<double>& parameters, 
  const std::vector<double>& temperature, 
  const std::vector<double>& pressure,
  std::vector<std::vector<double>>& number_densities, 
  std::vector<double>& mean_molecular_weight)
{
  //first, set the total number density (in case we're the first chemistry model)
  for (size_t i=0; i<number_densities.size(); ++i)
    number_densities[i][_TOTAL] = pressure[i] * 1.e6 / constants::boltzmann_k / temperature[i];

  const double mixing_ratio_1 = parameters[0];
  const double mixing_ratio_2 = parameters[1];
  const double pressure_transition = parameters[2];

  bool neglect_model = false;

  unsigned int transition_index = 0;

  for (size_t i=0; i<pressure.size(); ++i)
  {
    if (pressure[i] < pressure_transition)
    {
      transition_index = i;
      break;
    }
  }

  for (size_t i=0; i<transition_index; ++i)
    number_densities[i][species.front()] = number_densities[i][_TOTAL] * mixing_ratio_1;
  
  for (size_t i=transition_index; i<pressure.size(); ++i)
    number_densities[i][species.front()] = number_densities[i][_TOTAL] * mixing_ratio_2;

  //calculate the mean molecular weight
  //note that we sum over *all* species, not just the ones that were included in this chemistry
  meanMolecularWeight(number_densities, mean_molecular_weight);

  bool mixing_ratios_ok = checkMixingRatios(number_densities);

  if (!mixing_ratios_ok)
    neglect_model = true;

  return neglect_model;
}



}

