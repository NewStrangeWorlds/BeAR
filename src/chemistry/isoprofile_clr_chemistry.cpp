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

#include "isoprofile_clr_chemistry.h"

#include "chem_species.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"


namespace helios {


IsoprofileCLRChemistry::IsoprofileCLRChemistry(const std::vector<std::string>& chemical_species)
{ 
  std::cout << "- Chemistry model: " << "isoprofiles, centred-log-ratio prior" << "\n";
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

  nb_parameters = chemical_species.size() - 1;
}




bool IsoprofileCLRChemistry::calcChemicalComposition(
  const std::vector<double>& parameters,
  const std::vector<double>& temperature,
  const std::vector<double>& pressure,
  std::vector<std::vector<double>>& number_densities,
  std::vector<double>& mean_molecular_weight)
{
  bool neglect_model = false;

  std::vector<double> mixing_ratios(species.size(), 0.0);

  neglect_model = transformPriors(parameters, mixing_ratios);


  for (size_t i=0; i<number_densities.size(); ++i)
  { 
    number_densities[i][_TOTAL] = pressure[i] * 1.e6 / constants::boltzmann_k / temperature[i];

    for (size_t j=0; j<species.size(); ++j)
      number_densities[i][species[j]] = number_densities[i][_TOTAL] * mixing_ratios[j];
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


//transform priors to mixing ratios
//returns bool if the model is to be neglected
bool IsoprofileCLRChemistry::transformPriors(
  const std::vector<double>& priors,
  std::vector<double>& mixing_ratios)
{
  const size_t nb_species = species.size();

  const double ln_g_min = 1.0/nb_species 
    * (std::log(min_number_density) 
    + (nb_species - 1.0) * std::log((1.0 - min_number_density)/(nb_species - 1.0)));

  const double xi_min = std::log(min_number_density) - ln_g_min;

  const double ln_g_max = 1.0/nb_species 
    * (std::log(1.0 - (nb_species-1.0)*min_number_density) 
    + (nb_species - 1.0)*std::log(min_number_density) );
  
  const double xi_max = std::log(1.0 - (nb_species-1.0)*min_number_density) - ln_g_max;
  
  std::vector<double> transformed_priors(nb_species, 0.0);

  for (size_t i=0; i<nb_species-1; ++i)
  {
    transformed_priors[i] = priors[i] * (xi_max - xi_min) + xi_min;
    transformed_priors[nb_species-1] -= transformed_priors[i];
  }

  if (transformed_priors[nb_species-1] > xi_max)
    return true;


  double normalisation = 0;

  for (auto & p : transformed_priors)
    normalisation += std::exp(p);


  for (size_t i=0; i<nb_species; ++i)
    mixing_ratios[i] = std::exp(transformed_priors[i])/normalisation;


  if (*std::min_element(mixing_ratios.begin(), mixing_ratios.end()) < min_number_density)
    return true;


  return false;
}



}

