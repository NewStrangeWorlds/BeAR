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
#include <cmath>
#include <vector>
#include <string>

#include "iso_ratio_chemistry.h"

#include "chem_species.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"


namespace bear {


IsoprofileRatioChemistry::IsoprofileRatioChemistry(
  const std::vector<std::string>& chemical_species)
{
  if (chemical_species.size() < 2)
  {
    std::string error_message =
      "iso_ratio chemistry requires at least two species "
      "(one parent + one isotopologue).\n";
    throw InvalidInput(
      std::string("IsoprofileRatioChemistry::IsoprofileRatioChemistry"),
      error_message);
  }

  std::cout << "- Chemistry model: isoprofile with isotopologue ratio\n";
  std::cout << "  - Species: ";
  for (auto& s : chemical_species) std::cout << s << "  ";
  std::cout << "\n";
  std::cout << "  - Terrestrial standard ratio (parent/iso): " << standard_ratio << "\n";

  for (auto& name : chemical_species)
  {
    bool found = false;
    for (size_t j = 0; j < constants::species_data.size(); ++j)
    {
      if (constants::species_data[j].symbol == name)
      {
        species.push_back(constants::species_data[j].id);
        found = true;
        break;
      }
    }
    if (!found)
    {
      std::string error_message =
        "Chemical species " + name +
        " not found in the list of species in chem_species.h\n";
      throw InvalidInput(
        std::string("IsoprofileRatioChemistry::IsoprofileRatioChemistry"),
        error_message);
    }
  }

  //parameter_names is the source of truth for the parameter count:
  //N species -> N parameters, the first N-1 being direct isoprofile mixing
  //ratios and the last a log10 offset from the terrestrial isotopologue ratio
  for (size_t i = 0; i < chemical_species.size(); ++i)
  {
    std::string name = chemical_species[i];
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);

    if (i + 1 < chemical_species.size())
      parameter_names.push_back("chem_" + name + "_mr");
    else
      parameter_names.push_back("chem_" + name + "_ratio");
  }
}



bool IsoprofileRatioChemistry::calcChemicalComposition(
  const std::vector<double>& parameters,
  const std::vector<double>& temperature,
  const std::vector<double>& pressure,
  std::vector<std::vector<double>>& number_densities,
  std::vector<double>& mean_molecular_weight)
{
  const size_t nb_direct = species.size() - 1;   // all but the last species
  const size_t parent_idx = nb_direct - 1;        // second-to-last species

  // Literature convention (Line et al. 2021; Smith et al. 2023):
  //   [13CO/12CO] = log(13CO/12CO)_planet - log(13CO/12CO)_Earth,
  //   (13CO/12CO)_Earth = 1/89  (Meibom et al. 2007)
  // so the retrieved parameter is a log10 OFFSET and
  //   13CO/12CO = (1/89) * 10^[13CO/12CO]
  // i.e. a POSITIVE value means more 13CO than terrestrial. An earlier version
  // divided by 10^p instead of multiplying, which silently flipped the sign:
  // [13CO/12CO] = +0.5 gave 1/281 where the literature gives 1/28.
  const double iso_vmr =
    parameters[parent_idx] / standard_ratio * std::pow(10.0, parameters.back());

  bool neglect_model = false;

  for (size_t i = 0; i < number_densities.size(); ++i)
  {
    number_densities[i][_TOTAL] =
      pressure[i] * 1.e6 / constants::boltzmann_k / temperature[i];

    // Direct isoprofile species (H2O, CO, ...)
    for (size_t j = 0; j < nb_direct; ++j)
      number_densities[i][species[j]] = number_densities[i][_TOTAL] * parameters[j];

    // Derived isotopologue (13CO)
    number_densities[i][species.back()] = number_densities[i][_TOTAL] * iso_vmr;
  }

  meanMolecularWeight(number_densities, mean_molecular_weight);

  if (!checkMixingRatios(number_densities))
    neglect_model = true;

  return neglect_model;
}


}
