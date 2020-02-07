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


#include "piecewise_poly_temperature.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"

#include <algorithm>
#include <vector>
#include <cmath>


namespace helios {


PiecewisePolynomialTemperature::PiecewisePolynomialTemperature(const size_t nb_elements_in, const size_t polynomial_degree_in, const double atmos_boundaries [2])
 : temperature_profile(nb_elements_in, polynomial_degree_in, atmos_boundaries)
 , nb_elements{nb_elements_in}, polynomial_degree{polynomial_degree_in}
{

  nb_parameters = nb_elements*polynomial_degree + 1; //total number of temperature parameters

}



//calculate the temperature by a piecewise polynomial 
//the temperature will be evaluated on all points given by the pressure vector
bool PiecewisePolynomialTemperature::calcProfile(const std::vector<double>& parameters, const std::vector<double>& pressure, std::vector<double>& temperature)
{
  temperature.assign(pressure.size(), 0);

  //the temperature values at the degrees of freedom
  std::vector<double> temperature_dof(nb_parameters, 0.0);

  temperature_dof[0] = parameters[0];

  for (size_t i=1; i<nb_parameters; ++i)
    temperature_dof[i] = temperature_dof[i-1] * parameters[i];

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

