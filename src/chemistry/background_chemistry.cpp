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

#include "background_chemistry.h"

#include "chem_species.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"


namespace bear {


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

  if (chemical_species == "HHeEquilibrium")
  {
    hhe_eq_background = true;

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
  bool neglect_model = false;
  
  //first, set the total number density (in case we're the first chemistry model)
  for (size_t i=0; i<number_densities.size(); ++i)
    number_densities[i][_TOTAL] = pressure[i] * 1.e6 / constants::boltzmann_k / temperature[i];

  for (size_t i=0; i<number_densities.size(); ++i)
  { 
    //the mixing ratio of the background gas
    double mixing_ratio_background = mixingRatioBackground(number_densities[i]);
    
    //if we have a negative mixing ratio of the background gas
    //(i.e. the sum of all other species is larger than 1)
    //neglect this model
    if (mixing_ratio_background < 0)
    {
      neglect_model = true;
      mixing_ratio_background = 0;
    }
  
    if (h2he_background)
    {
      backgroundH2He(mixing_ratio_background, number_densities[i]);
    }
    else if (hhe_eq_background)
    {
      equilibriumHHe(pressure[i], temperature[i], mixing_ratio_background, number_densities[i]);
    }
    else
      number_densities[i][species.front()] = number_densities[i][_TOTAL] * mixing_ratio_background;
  }

  //calculate the mean molecular weight
  //note that we sum over *all* species, not just the ones that were included in this chemistry
  meanMolecularWeight(number_densities, mean_molecular_weight);
  
  return neglect_model;
}


double BackgroundChemistry::mixingRatioBackground(
  const std::vector<double>& number_densities)
{
  //the mixing ratio of the background gas
  double mixing_ratio_background = 1.0;

  for (auto & j : constants::species_data)
  {
    if (j.id != _TOTAL)
    {
      if (h2he_background)
      {
        if (j.id != _H2 && j.id != _He)
          mixing_ratio_background -= number_densities[j.id]/number_densities[_TOTAL];
      }
      else if (hhe_eq_background)
      {
        if (j.id != _H2 && j.id != _He && j.id != _H)
          mixing_ratio_background -= number_densities[j.id]/number_densities[_TOTAL];
      }
      else if (j.id != species.front())
        mixing_ratio_background -= number_densities[j.id]/number_densities[_TOTAL];
    }
  }

  return mixing_ratio_background;
}



void BackgroundChemistry::backgroundH2He(
  const double mixing_ratio_background, 
  std::vector<double>& number_densities)
{
  const double solar_h2 = 0.5;
  const double solar_he = 0.085114;
  const double solar_h2_he = solar_h2 + solar_he;
  //const double solar_h_he = solar_he + 1.0;
  //const double epsilon_h = 1.0 / (solar_h_he);
  
  //const double mixing_ratio_h = number_densities[i][_H]/number_densities[i][_TOTAL];
  //const double mixing_ratio_h2 = epsilon_h / (2.0 - epsilon_h) * (mixing_ratio_background - mixing_ratio_h);
  //const double mixing_ratio_he = mixing_ratio_background - mixing_ratio_h - mixing_ratio_h2;
  const double mixing_ratio_h2 = mixing_ratio_background * solar_h2 / solar_h2_he;
  const double mixing_ratio_he = mixing_ratio_background * solar_he / solar_h2_he;

  number_densities[_H2] = number_densities[_TOTAL] * mixing_ratio_h2;
  number_densities[_He] = number_densities[_TOTAL] * mixing_ratio_he;
}



//Computes the equilibrium abundances of H, H2, and He for a given pressure, temperature, 
//It distributes the background gas (with a given mixing ratio) between H, H2, and He according to chemical equilibrium. 
//The equilibrium constant for the reaction H2 -> 2H is taken from FastChem
void BackgroundChemistry::equilibriumHHe(
  const double pressure, 
  const double temperature,
  const double mixing_ratio_background, 
  std::vector<double>& number_densities)
{ //solar elemental abundances
  constexpr double X = 0.7060;            // H  mass fraction
  constexpr double Y = 0.2753;            // He mass fraction
  constexpr double HE_H = (Y / 4.0) / X;  // He/H_nuc ≈ 0.097486
  
  const double pressure_background = pressure * mixing_ratio_background;

  //lambda function for the equilibrium constant of the reaction H2 <-> 2H
   auto Kp_bar = [](double t) {
    constexpr double a1 =  5.1909637142380554e+04; //FastChem fit 
    constexpr double a2 = -1.8011701211306956e+00;
    constexpr double a3 =  8.7224583233705744e-02;
    constexpr double a4 =  2.5613890164973008e-04;
    constexpr double a5 = -5.3540255367406060e-09;

    const double ln_Kbar_form = a1/t + a2*std::log(t) + a3 + a4*t + a5*t*t;
    return std::exp(-ln_Kbar_form);
  };

  const double f  = HE_H;
  const double Kp = Kp_bar(temperature);
 
  const double A =  2.0*pressure_background + 0.5*Kp;
  const double B =  Kp*(f + 0.5);
  const double C = -Kp*(f + 1.0);
 
  double disc = B*B - 4.0*A*C;
  if (disc < 0.0) disc = 0.0;
 
  double alpha = (-B + std::sqrt(disc)) / (2.0*A);
  alpha = std::max(0.0, std::min(1.0, alpha));
 
  const double denom = f + (1.0 + alpha)/2.0;
  
  const double h_mixing_ratio = alpha/denom;
  const double h2_mixing_ratio = (1.0-alpha)/2.0/denom;
  const double he_mixing_ratio = f/denom;

  number_densities[_H] = number_densities[_TOTAL] * mixing_ratio_background * h_mixing_ratio;
  number_densities[_H2] = number_densities[_TOTAL] * mixing_ratio_background * h2_mixing_ratio;
  number_densities[_He] = number_densities[_TOTAL] * mixing_ratio_background * he_mixing_ratio;
}



}

