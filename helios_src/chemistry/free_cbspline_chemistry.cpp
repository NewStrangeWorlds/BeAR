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

#include "free_cbspline_chemistry.h"

#include "chem_species.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"

#include "../../_deps/boost_math-src/include/boost/math/interpolators/cardinal_cubic_b_spline.hpp"


namespace helios {


FreeCBSplineChemistry::FreeCBSplineChemistry(
  const std::string& chemical_species,
  const size_t nb_control_points_)
    : nb_control_points{nb_control_points_}
{ 
  std::cout << "- Chemistry model: " << "free cubic b-spline chemistry" << "\n";
  std::cout << "  - Species for this model: " << chemical_species << "\n";
  std::cout << "  - number of control points: " << nb_control_points << "\n";
  std::cout << "\n";


  if (nb_control_points < 5)
  {
    std::string error_message = "The free chemistry with cubic b splines requires at least 5 control points!\n";
    throw InvalidInput(std::string ("FreeCBSplineChemistry::FreeCBSplineChemistry"), error_message);
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
    throw InvalidInput(std::string ("FreeCBSplineChemistry::FreeCBSplineChemistry"), error_message);
  }


  nb_parameters = nb_control_points; //total number of chemistry parameters
}



//the parameters should contain the mixing ratios at the dof
bool FreeCBSplineChemistry::calcChemicalComposition(
  const std::vector<double>& parameters, 
  const std::vector<double>& temperature, 
  const std::vector<double>& pressure,
  std::vector<std::vector<double>>& number_densities, 
  std::vector<double>& mean_molecular_weight)
{
  //first, set the total number density (in case we're the first chemistry model)
  for (size_t i=0; i<number_densities.size(); ++i)
    number_densities[i][_TOTAL] = pressure[i] * 1.e6 / constants::boltzmann_k / temperature[i];


  if (parameters.size() != nb_control_points)
    std::cout << "The number of free parameters is not equal to the number of control points!\n";


  double control_points_step = (std::log10(pressure[0]) - std::log10(pressure.back())) / (nb_control_points - 1.0);


  std::vector<double> mixing_ratio_control_point(nb_control_points, 0.0);

  for (size_t i=0; i<nb_parameters; ++i)
    mixing_ratio_control_point[i] = std::log10(parameters[i]);

  std::reverse(mixing_ratio_control_point.begin(), mixing_ratio_control_point.end());

  boost::math::interpolators::cardinal_cubic_b_spline<double> mixing_ratios (
    mixing_ratio_control_point.data(), 
    mixing_ratio_control_point.size(), 
    std::log10(pressure.back()), 
    control_points_step);


  bool neglect_model = false;

  //calculate the number density at each level from the cubic spline
  //also, check for inconsistencies
  for (size_t i=0; i<pressure.size(); ++i)
  {
    double mixing_ratio = std::pow(10, mixing_ratios(std::log10(pressure[i])));

    if (mixing_ratio > 1)
    {
      mixing_ratio = 1;
      neglect_model = true;
    }

    number_densities[i][species.front()] = number_densities[i][_TOTAL] * mixing_ratio;

    //adjust the mean molecular weight
    mean_molecular_weight[i] += 
      number_densities[i][species.front()]/number_densities[i][_TOTAL] 
      * constants::species_data[species.front()].molecular_weight;
  }


  return neglect_model;
}



}

