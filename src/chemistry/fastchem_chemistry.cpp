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


#include "fastchem_chemistry.h"

#include "chem_species.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"
#include "../../_deps/fastchem-src/fastchem_src/fastchem.h"


#include <algorithm>
#include <unordered_map>
#include <vector>
#include <omp.h>
#include <math.h>
#include <iostream>
#include <string>


namespace bear {


FastChemChemistry::FastChemChemistry(
  const std::string& fastchen_parameter_file,
  const size_t nb_openmp_proc,
  const std::vector<std::string>& ratio_specs)
  : fastchem(fastchen_parameter_file, 1)
  , nb_processes{nb_openmp_proc}
{
  std::cout << "- Chemistry model: " << "equilibrium/FastChem" << "\n";
  std::cout << "  - Parameter file: " << fastchen_parameter_file << "\n";

  reference_element_abundances = fastchem.getElementAbundances();
  //std::cout << fastchem.getGasSpeciesIndex("H2O1888") << "\n"; exit(0);
  // Build reverse lookup: FastChem species symbol -> FastChem index.
  // We iterate by index rather than calling getGasSpeciesIndex (which uses
  // find_if over a pointer vector and triggers an optimizer crash at -O3).
  const unsigned int nb_fc_species = fastchem.getGasSpeciesNumber();
  std::unordered_map<std::string, size_t> fc_index_map;
  for (unsigned int j = 0; j < nb_fc_species; ++j)
    fc_index_map[fastchem.getGasSpeciesSymbol(j)] = j;

  fastchem_species_indices.assign(constants::species_data.size(), fastchem::FASTCHEM_UNKNOWN_SPECIES);
  for (size_t i = 0; i < constants::species_data.size(); ++i)
  {
    auto it = fc_index_map.find(constants::species_data[i].fastchem_symbol);
    if (it != fc_index_map.end())
      fastchem_species_indices[i] = it->second;
  }

  if (fastchem_species_indices[_H] == fastchem::FASTCHEM_UNKNOWN_SPECIES
      || fastchem_species_indices[_O] == fastchem::FASTCHEM_UNKNOWN_SPECIES
      || fastchem_species_indices[_C] == fastchem::FASTCHEM_UNKNOWN_SPECIES)
  {
    std::string error_message = "Critical elements (H, C, or O) not found in FastChem\n";
    throw InvalidInput(std::string ("FastChemChemistry::FastChemChemistry"), error_message);
  }

  for (const auto& spec : ratio_specs)
  {
    const size_t slash = spec.find('/');
    if (slash == std::string::npos || slash == 0 || slash == spec.size() - 1)
    {
      std::string error_message = "Invalid element ratio specification '" + spec + "'; expected format X/Y\n";
      throw InvalidInput(std::string ("FastChemChemistry::FastChemChemistry"), error_message);
    }

    const std::string num_symbol = spec.substr(0, slash);
    const std::string den_symbol = spec.substr(slash + 1);

    const unsigned int num_idx = fastchem.getElementIndex(num_symbol);
    const unsigned int den_idx = fastchem.getElementIndex(den_symbol);

    if (num_idx == fastchem::FASTCHEM_UNKNOWN_SPECIES)
    {
      std::string error_message = "Element '" + num_symbol + "' from ratio '" + spec + "' not found in FastChem\n";
      throw InvalidInput(std::string ("FastChemChemistry::FastChemChemistry"), error_message);
    }
    if (den_idx == fastchem::FASTCHEM_UNKNOWN_SPECIES)
    {
      std::string error_message = "Element '" + den_symbol + "' from ratio '" + spec + "' not found in FastChem\n";
      throw InvalidInput(std::string ("FastChemChemistry::FastChemChemistry"), error_message);
    }

    const unsigned int h_idx = fastchem.getElementIndex("H");
    const double ref_ratio = (den_idx == h_idx)
        ? reference_element_abundances[num_idx] / reference_element_abundances[den_idx]
        : 1.0;
    element_ratios.push_back({num_idx, den_idx, ref_ratio, spec});
    std::cout << "  - Element ratio: " << spec << "\n";
  }

  std::cout << "\n";

  nb_parameters = 1 + element_ratios.size();
}




bool FastChemChemistry::calcChemicalComposition(
  const std::vector<double>& parameters,
  const std::vector<double>& temperature,
  const std::vector<double>& pressure,
  std::vector<std::vector<double>>& number_densities,
  std::vector<double>& mean_molecular_weight)
{
  const double metallicity_factor = parameters[0];

  std::vector<double> element_abundances = reference_element_abundances;

  for (size_t i=0; i<element_abundances.size(); ++i)
    if (i != fastchem_species_indices[_H] && (i != fastchem_species_indices[_He] && fastchem_species_indices[_He] != fastchem::FASTCHEM_UNKNOWN_SPECIES) )
      element_abundances[i] *= metallicity_factor;

  for (size_t i=0; i<element_ratios.size(); ++i)
    element_abundances[element_ratios[i].numerator_idx] =
      element_abundances[element_ratios[i].denominator_idx]
      * parameters[1 + i]
      * element_ratios[i].reference_ratio;


  fastchem.setElementAbundances(element_abundances);


  //set up the input & output structures and run the chemistry
  fastchem::FastChemInput input;
  fastchem::FastChemOutput output;

  input.temperature = temperature;
  input.pressure = pressure;

  for (auto & i : input.temperature)
    if (i < 100) i = 100;


  bool neglect_model = false;


  size_t status = fastchem.calcDensities(input, output);

  if (status == fastchem::FASTCHEM_INITIALIZATION_FAILED)
  {
    std::string error_message = "FastChem initialisation failed!\n";
    throw InvalidInput(std::string ("FastChemChemistry::calcChemicalComposition"), error_message);
  }

  if (status != fastchem::FASTCHEM_SUCCESS)
    neglect_model = true;

  
  mean_molecular_weight = output.mean_molecular_weight;

  for (size_t j=0; j<temperature.size(); ++j)
    for (size_t i=0; i<constants::species_data.size(); ++i)
    {
      if (fastchem_species_indices[i] != fastchem::FASTCHEM_UNKNOWN_SPECIES)
        number_densities[j][i] = output.number_densities[j][fastchem_species_indices[i]];

      if (i == _TOTAL)
        number_densities[j][_TOTAL] = pressure[j] * 1.e6 / constants::boltzmann_k / temperature[j];
    }

  return neglect_model;
}



}

