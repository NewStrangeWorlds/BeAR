// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include "piecewise_poly_temperature.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"

#include <algorithm>
#include <vector>
#include <cmath>
#include <string>


namespace bear {


PiecewisePolynomialTemperature::PiecewisePolynomialTemperature(
  const size_t nb_elements_in,
  const size_t polynomial_degree_in,
  const std::vector<double>& atmos_boundaries,
  const Parametrisation parametrisation_)
 : temperature_profile(nb_elements_in, polynomial_degree_in, atmos_boundaries)
 , nb_elements{nb_elements_in}, polynomial_degree{polynomial_degree_in}
{
  parametrisation = parametrisation_;

  //parameter_names is the source of truth for the parameter count.
  const size_t nb_dof = nb_elements*polynomial_degree + 1;

  setControlPointNames(nb_dof);
}



//calculate the temperature by a piecewise polynomial 
//the temperature will be evaluated on all points given by the pressure vector
bool PiecewisePolynomialTemperature::calcProfile(
  const std::vector<double>& parameters,
  const double surface_gravity,
  const std::vector<double>& pressure,
  std::vector<double>& temperature)
{
  temperature.assign(pressure.size(), 0);

  //the temperature values at the degrees of freedom
  std::vector<double> temperature_dof = controlTemperatures(parameters);

  temperature_profile.setDOFvalues(temperature_dof);


  for (size_t i=0; i<pressure.size(); ++i)
  {
    double p_log = std::log10(pressure[i]);

    temperature[i] = temperature_profile.getValue(p_log);
  }

  
  //neglect models with too low temperatures
  bool neglect_model = false;
  
  for (auto & i : temperature)
    if (i < 50) {i = 50; neglect_model = true;}


  return neglect_model;
}


}

