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

#include <algorithm>
#include <vector>
#include <cmath>

#include "free_chemistry.h"

#include "chem_species.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"
#include "../additional/piecewise_poly.h"


namespace helios {


FreeChemistry::FreeChemistry(
  const std::string& chemical_species,
  const size_t nb_elements_in, 
  const size_t polynomial_degree_in, 
  const double atmos_boundaries [2])
    : mixing_ratios(nb_elements_in, polynomial_degree_in, atmos_boundaries)
    , nb_elements{nb_elements_in}, polynomial_degree{polynomial_degree_in}
{ 
  std::cout << "- Chemistry model: " << "free chemistry" << "\n";
  std::cout << "  - Species for this model: " << chemical_species << "\n";
  std::cout << "  - number of elements: " << nb_elements_in << "\n";
  std::cout << "  - polynomial degree: " << polynomial_degree_in << "\n";
  std::cout << "\n";


  if (nb_elements < 1 || polynomial_degree < 1)
  {
    std::string error_message = "Wrong number of elements or polynomial degree in free chemistry config\n";
    throw InvalidInput(std::string ("FreeChemistry::FreeChemistry"), error_message);
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
    throw InvalidInput(std::string ("FreeChemistry::FreeChemistry"), error_message);
  }
      

  //std::cout << "Free chemistry initialised with the following species: \n";  
  //std::cout << constants::species_data[species.front()].symbol << "\n";
  

  nb_parameters = nb_elements*polynomial_degree + 1; //total number of chemistry parameters
}



//the parameters should contain the mixing ratios at the dof
bool FreeChemistry::calcChemicalComposition(
  const std::vector<double>& parameters, 
  const std::vector<double>& temperature, 
  const std::vector<double>& pressure,
  std::vector<std::vector<double>>& number_densities, 
  std::vector<double>& mean_molecular_weight)
{
  //first, set the total number density (in case we're the first chemistry model)
  for (size_t i=0; i<number_densities.size(); ++i)
    number_densities[i][_TOTAL] = pressure[i] * 1.e6 / constants::boltzmann_k / temperature[i];


  std::vector<double> mixing_ratios_dof = parameters;
  
  for (auto & i : mixing_ratios_dof)
    i = std::log10(i);

  mixing_ratios.setDOFvalues(mixing_ratios_dof);


  bool neglect_model = false;

  //calculate the number density at each level from the piecewise polynomial
  //also, check for inconsistencies
  for (size_t i=0; i<pressure.size(); ++i)
  {
    const double p_log = std::log10(pressure[i]);
    double mixing_ratio = std::pow(10, mixing_ratios.getValue(p_log));

    if (mixing_ratio > 1)
    {
      mixing_ratio = 1;
      neglect_model = true;
    }

    number_densities[i][species.front()] = number_densities[i][_TOTAL] * mixing_ratio;
    
    //adjust the mean molecular weight
    mean_molecular_weight[i] += number_densities[i][species.front()]/number_densities[i][_TOTAL] * constants::species_data[species.front()].molecular_weight;
  }
  

  return neglect_model;
}



}

