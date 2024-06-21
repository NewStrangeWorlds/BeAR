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

#include "background_chemistry.h"

#include "chem_species.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"


namespace helios {


BackgroundChemistry::BackgroundChemistry(const std::string& chemical_species)
{ 
  std::cout << "- Chemistry model: " << "background" << "\n";
  std::cout << "  - Species for this model: " << chemical_species << "\n";
  
  nb_parameters = 0;


  if (chemical_species == "H2He")
  {
    h2he_background = true;

    return;
  }


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
    throw InvalidInput(std::string ("BackgroundChemistry::BackgroundChemistry"), error_message);
  }

}



bool BackgroundChemistry::calcChemicalComposition(
  const std::vector<double>& parameters,
  const std::vector<double>& temperature,
  const std::vector<double>& pressure,
  std::vector<std::vector<double>>& number_densities,
  std::vector<double>& mean_molecular_weight)
{
  const double solar_h2_he = solar_h2 + solar_he;
  //const double solar_h_he = solar_he + 1.0;
  //const double epsilon_h = 1.0 / (solar_h_he);

  bool neglect_model = false;

  //first, set the total number density (in case we're the first chemistry model)
  for (size_t i=0; i<number_densities.size(); ++i)
    number_densities[i][_TOTAL] = pressure[i] * 1.e6 / constants::boltzmann_k / temperature[i];


  for (size_t i=0; i<number_densities.size(); ++i)
  { 
    //the mixing ratio of the background gas
    double mixing_ratio_background = 1.0;

    for (auto & j : constants::species_data)
      if (j.id != _TOTAL)
      {
        if (h2he_background)
        {
          if (j.id != _H2 && j.id != _He)
             mixing_ratio_background -= number_densities[i][j.id]/number_densities[i][_TOTAL];
        }
        else if (j.id != species.front())
          mixing_ratio_background -= number_densities[i][j.id]/number_densities[i][_TOTAL];
      }

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
    if (h2he_background)
    {
      const double mixing_ratio_h2 = mixing_ratio_background * solar_h2 / solar_h2_he;
      const double mixing_ratio_he = mixing_ratio_background * solar_he / solar_h2_he;

      number_densities[i][_H2] = number_densities[i][_TOTAL] * mixing_ratio_h2;
      number_densities[i][_He] = number_densities[i][_TOTAL] * mixing_ratio_he;
      continue;
    }
    
    number_densities[i][species.front()] = number_densities[i][_TOTAL] * mixing_ratio_background;
  }

  //calculate the mean molecular weight
  //note that we sum over *all* species, not just the ones that were included in this chemistry
  meanMolecularWeight(number_densities, mean_molecular_weight);

  return neglect_model;
}



}

