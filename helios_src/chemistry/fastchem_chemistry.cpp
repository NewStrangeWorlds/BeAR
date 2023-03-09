/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2022 Daniel Kitzmann
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


#include "fastchem_chemistry.h"

#include "chem_species.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"
#include "../../_deps/fastchem-src/fastchem_src/fastchem.h"


#include <algorithm>
#include <vector>
#include <omp.h>
#include <math.h> 
#include <iostream>


namespace helios {


FastChemChemistry::FastChemChemistry(
  const std::string& fastchen_parameter_file, const size_t nb_openmp_proc)
  : fastchem(fastchen_parameter_file, 1) 
  , nb_processes{nb_openmp_proc}
{
  std::cout << "- Chemistry model: " << "equilibrium/FastChem" << "\n";
  std::cout << "  - Parameter file: " << fastchen_parameter_file << "\n\n";
  

  reference_element_abundances = fastchem.getElementAbundances();

  
  fastchem_species_indices.assign(constants::species_data.size(), fastchem::FASTCHEM_UNKNOWN_SPECIES);

  for (size_t i=0; i<constants::species_data.size(); ++i)
    fastchem_species_indices[i] = fastchem.getSpeciesIndex(constants::species_data[i].fastchem_symbol);
    
  
  //check if C, O, and H are present in FastChem
  if (fastchem_species_indices[_H] == fastchem::FASTCHEM_UNKNOWN_SPECIES 
      || fastchem_species_indices[_O] == fastchem::FASTCHEM_UNKNOWN_SPECIES 
      || fastchem_species_indices[_C] == fastchem::FASTCHEM_UNKNOWN_SPECIES)
  {
    std::string error_message = "Critical elements (H, C, or O) not found in FastChem\n";
    throw InvalidInput(std::string ("FastChemChemistry::FastChemChemistry"), error_message);
  }
  

  nb_parameters = 2;
}




bool FastChemChemistry::calcChemicalComposition(
  const std::vector<double>& parameters,
  const std::vector<double>& temperature,
  const std::vector<double>& pressure,
  std::vector<std::vector<double>>& number_densities,
  std::vector<double>& mean_molecular_weight)
{
  const double metallicity_factor = parameters[0];
  const double co_ratio = parameters[1];

  std::vector<double> element_abundances = reference_element_abundances;
  
  for (size_t i=0; i<element_abundances.size(); ++i)
    if (i != fastchem_species_indices[_H] && (i != fastchem_species_indices[_He] && fastchem_species_indices[_He] != fastchem::FASTCHEM_UNKNOWN_SPECIES) )
      element_abundances[i] *= metallicity_factor;


  element_abundances[fastchem_species_indices[_O]] = element_abundances[fastchem_species_indices[_C]] / co_ratio;


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

