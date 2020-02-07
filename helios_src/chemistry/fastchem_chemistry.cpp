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


#include "fastchem_chemistry.h"

#include "chem_species.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"
#include "../../fastchem_src/fastchem.h"


#include <algorithm>
#include <vector>
#include <omp.h>



namespace helios {


FastChemChemistry::FastChemChemistry(const std::string& fastchen_parameter_file, const size_t nb_openmp_proc)
  : fastchem(fastchen_parameter_file, 1) 
  , nb_processes{nb_openmp_proc}
{
  
  reference_element_abundances = fastchem.getElementAbundance();

  
  fastchem_species_indices.assign(constants::species_data.size(), fastchem::FASTCHEM_UNKNOWN_SPECIES);

  for (size_t i=0; i<constants::species_data.size(); ++i)
    fastchem_species_indices[i] = fastchem.getSpeciesIndex(constants::species_data[i].fastchem_symbol);
    
  
  //check if C, O, and H are present in FastChem
  if (fastchem_species_indices[_H] == fastchem::FASTCHEM_UNKNOWN_SPECIES 
      || fastchem_species_indices[_O] == fastchem::FASTCHEM_UNKNOWN_SPECIES 
      || fastchem_species_indices[_C] == fastchem::FASTCHEM_UNKNOWN_SPECIES)
  {
    std::string error_message = "Critical elements (H, C, or O) not found in FastChem\n";
    throw ExceptionInvalidInput(std::string ("FastChemChemistry::FastChemChemistry"), error_message);
  }
  

  nb_parameters = 2;
}




bool FastChemChemistry::calcChemicalComposition(const std::vector<double>& parameters, const std::vector<double>& temperature, const std::vector<double>& pressure,
                                                  std::vector<std::vector<double>>& number_densities, std::vector<double>& mean_molecular_weight)
{
  const double metallicity_factor = parameters[0];
  const double co_ratio = parameters[1];

  std::vector<double> element_abundances = reference_element_abundances;
  
  for (size_t i=0; i<element_abundances.size(); ++i)
    if (i != fastchem_species_indices[_H] && (i != fastchem_species_indices[_He] && fastchem_species_indices[_He] != fastchem::FASTCHEM_UNKNOWN_SPECIES) )
      element_abundances[i] *= metallicity_factor;


  element_abundances[fastchem_species_indices[_O]] = element_abundances[fastchem_species_indices[_C]] / co_ratio;
  
  
  fastchem.setElementAbundance(element_abundances);

  std::vector< fastchem::FastChem<long double>  > fastchems(nb_processes, fastchem);


  std::vector<double> h_densities(temperature.size(), 0.0);
  std::vector<std::vector<double>> fastchem_results;
  fastchem_results.resize(temperature.size());

  
  bool neglect_model = false;

  #pragma omp parallel for schedule(dynamic, 1)
  for (unsigned int i=0; i<temperature.size(); i++)
  {
    double local_temperature = temperature[i];
    const double local_pressure = pressure[i] * 1e6;

    if (local_temperature < 200) local_temperature = 200;

   
    size_t status = fastchems[omp_get_thread_num()].calcDensities(local_temperature, local_pressure, fastchem_results[i], h_densities[i], mean_molecular_weight[i]);

    if (status == fastchem::FASTCHEM_INITIALIZATION_FAILED)
    {
      std::string error_message = "FastChem initialisation failed!\n";
      throw ExceptionInvalidInput(std::string ("FastChemChemistry::calcChemicalComposition"), error_message);
    }


    if (status != fastchem::FASTCHEM_SUCCESS)
      neglect_model = true;
  }

  
  for (size_t j=0; j<temperature.size(); ++j)
    for (size_t i=0; i<constants::species_data.size(); ++i)
      if (fastchem_species_indices[i] != fastchem::FASTCHEM_UNKNOWN_SPECIES)
        number_densities[j][i] = fastchem_results[j][fastchem_species_indices[i]];
  

  return neglect_model;
}



}

