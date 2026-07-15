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
#include <string>

#include "free_pchip_chemistry.h"

#include "chem_species.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"

#include "../../_deps/boost_math-src/include/boost/math/interpolators/pchip.hpp"


namespace bear {


FreePchipChemistry::FreePchipChemistry(
  const std::string& chemical_species,
  const size_t nb_control_points_)
    : nb_control_points{nb_control_points_}
{
  std::cout << "- Chemistry model: " << "free PCHIP chemistry" << "\n";
  std::cout << "  - Species for this model: " << chemical_species << "\n";
  std::cout << "  - number of control points: " << nb_control_points << "\n";
  std::cout << "\n";

  if (nb_control_points < 4)
  {
    std::string error_message = "The free PCHIP chemistry requires at least 4 control points!\n";
    throw InvalidInput(std::string ("FreePchipChemistry::FreePchipChemistry"), error_message);
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
    throw InvalidInput(std::string ("FreePchipChemistry::FreePchipChemistry"), error_message);
  }

  //parameter_names is the source of truth for the parameter count:
  //one control point per parameter. Each parameter is the ABSOLUTE (log)
  //mixing ratio at a control point (control_point[i] = log10(parameters[i])),
  //hence the _t suffix.
  std::string species_name = chemical_species;
  std::transform(species_name.begin(), species_name.end(), species_name.begin(), ::tolower);

  for (size_t i=0; i<nb_control_points; ++i)
    parameter_names.push_back("mr_" + species_name + "_t" + std::to_string(i));
}



bool FreePchipChemistry::calcChemicalComposition(
  const std::vector<double>& parameters,
  const std::vector<double>& temperature,
  const std::vector<double>& pressure,
  std::vector<std::vector<double>>& number_densities,
  std::vector<double>& mean_molecular_weight)
{
  for (size_t i=0; i<number_densities.size(); ++i)
    number_densities[i][_TOTAL] = pressure[i] * 1.e6 / constants::boltzmann_k / temperature[i];

  if (parameters.size() != nb_control_points)
    std::cout << "The number of free parameters is not equal to the number of control points!\n";

  double control_points_step = (std::log10(pressure[0]) - std::log10(pressure.back())) / (nb_control_points - 1.0);

  std::vector<double> mixing_ratio_control_point(nb_control_points, 0.0);

  for (size_t i=0; i<nb_control_points; ++i)
    mixing_ratio_control_point[i] = std::log10(parameters[i]);

  std::reverse(mixing_ratio_control_point.begin(), mixing_ratio_control_point.end());

  std::vector<double> x_knots(nb_control_points);
  for (size_t i = 0; i < nb_control_points; ++i)
    x_knots[i] = std::log10(pressure.back()) + i * control_points_step;

  auto mixing_ratios = boost::math::interpolators::pchip<std::vector<double>>(
    std::move(x_knots), std::move(mixing_ratio_control_point));

  bool neglect_model = false;

  for (size_t i=0; i<pressure.size(); ++i)
  {
    double mixing_ratio = std::pow(10, mixing_ratios(std::log10(pressure[i])));

    if (mixing_ratio > 1)
    {
      mixing_ratio = 1;
      neglect_model = true;
    }

    number_densities[i][species.front()] = number_densities[i][_TOTAL] * mixing_ratio;
  }

  meanMolecularWeight(number_densities, mean_molecular_weight);

  bool mixing_ratios_ok = checkMixingRatios(number_densities);

  if (!mixing_ratios_ok)
    neglect_model = true;

  return neglect_model;
}


}
