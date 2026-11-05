// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <algorithm>
#include <vector>

#include "fixed_chemistry.h"

#include "chem_species.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"


namespace bear {


FixedChemistry::FixedChemistry(const std::vector<std::string>& species_symbol)
{

  findSpecies(species_symbol);

}


void FixedChemistry::findSpecies(
  const std::vector<std::string>& species_symbol)
{
  species.clear();

  for (auto & s : species_symbol)
  {
    bool species_found = false;

    for (size_t j=0; j<constants::species_data.size(); ++j)
    {
      if (constants::species_data[j].symbol == s || constants::species_data[j].fastchem_symbol == s)
      {
        species.push_back(constants::species_data[j].id); 
        species_found = true;
        break;
      }
    }

    if (!species_found)
    {
      std::string error_message = "Chemical species " + s + " not found in the list of species in chem_species.h \n";
      throw InvalidInput(std::string ("FixedChemistry::findSpecies"), error_message);
    }
  }
}


void FixedChemistry::setChemicalComposition(
  const std::vector<double>& pressure,
  const std::vector<double>& temperature,
  const std::vector<std::vector<double>>& mixing_ratios,
  std::vector<std::vector<double>>& number_densities,
  std::vector<double>& mean_molecular_weight)
{
  for (size_t i=0; i<number_densities.size(); ++i)
  { 
    number_densities[i][_TOTAL] = pressure[i] * 1.e6 / constants::boltzmann_k / temperature[i];

    for (size_t j=0; j<species.size(); ++j)
      number_densities[i][species[j]] = number_densities[i][_TOTAL] * mixing_ratios[i][j];
  }

  meanMolecularWeight(number_densities, mean_molecular_weight);
  
  bool mixing_ratios_ok = checkMixingRatios(number_densities);

  if (!mixing_ratios_ok)
    std::cout << "Warning! Sum of all mixing ratios is larger than 1.0! \n";
}


}