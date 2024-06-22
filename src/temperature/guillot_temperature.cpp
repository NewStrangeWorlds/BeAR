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


#include "guillot_temperature.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"

#include <algorithm>
#include <vector>
#include <cmath>


namespace bear {


GuillotTemperature::GuillotTemperature(const std::string profile_type)
{
  nb_parameters = 5;

  std::cout << "\n- Temperature profile: Guillot profile\n\n";
  
  if (profile_type == "beam")
  {
    profile = 0;
  }
  else if (profile_type == "isotropic")
  {
    profile = 1;
  }
  else
  {
    std::string error_message = "Profile type " + profile_type + " unknown!";
    throw InvalidInput(std::string ("GuillotTemperature::GuillotTemperature"), error_message);
  }

}


//calculate the temperature based on the analytical Milne solution
bool GuillotTemperature::calcProfile(
  const std::vector<double>& parameters, 
  const double surface_gravity,
  const std::vector<double>& pressure, 
  std::vector<double>& temperature)
{
  temperature.assign(pressure.size(), 0);

  const double kappa_ir = parameters[0];
  const double temperarure_irr = parameters[1];
  const double temperarure_int = parameters[2];
  const double gamma = parameters[3];
  const double profile_parameter = parameters[4];


  std::vector<double> tau_ir(pressure.size(), 0);

  for (size_t i=0; i<temperature.size(); ++i)
    tau_ir[i] = kappa_ir * pressure[i] * 1e6 / surface_gravity;


  if (profile == 0)
    profileBeamSource(
      tau_ir,
      temperarure_irr,
      temperarure_int,
      profile_parameter,
      gamma,
      temperature);

  if (profile == 1)
    profileIsotropicSource(
      tau_ir,
      temperarure_irr,
      temperarure_int,
      profile_parameter,
      gamma,
      temperature);

  //neglect models with too low temperatures
  bool neglect_model = false;
  
  for (auto & i : temperature)
    if (i < 50) {i = 50; neglect_model = true;}


  return neglect_model;
}


void GuillotTemperature::profileBeamSource(
      const std::vector<double>& optical_depth,
      const double temperature_irr,
      const double temperature_int,
      const double mu,
      const double gamma,
      std::vector<double>& temperature)
{
  for (size_t i=0; i<temperature.size(); ++i)
  {
    temperature[i] = 3.0/4.0 * std::pow(temperature_int, 4) * (2.0/3.0 + optical_depth[i])
                   + 3.0/4.0 * std::pow(temperature_irr, 4) * mu 
                   * (2.0/3.0 + mu/gamma + (gamma/(3.0 * mu) - mu/gamma) 
                      * std::exp(-gamma*optical_depth[i]/mu));

    temperature[i] = std::pow(temperature[i], 0.25);
  }

}


void GuillotTemperature::profileIsotropicSource(
      const std::vector<double>& optical_depth,
      const double temperature_irr,
      const double temperature_int,
      const double flux_distribution,
      const double gamma,
      std::vector<double>& temperature)
{
  for (size_t i=0; i<temperature.size(); ++i)
  {
    temperature[i] = 3.0/4.0 * std::pow(temperature_int, 4) * (2.0/3.0 + optical_depth[i])
                   + 3.0/4.0 * std::pow(temperature_irr, 4) * flux_distribution 
                   * (2.0/3.0 + 1.0/gamma/std::sqrt(3.0) + (gamma/std::sqrt(3.0) - 1.0/gamma/std::sqrt(3.0)) 
                     * std::exp(-gamma*optical_depth[i]*std::sqrt(3.0)));

    temperature[i] = std::pow(temperature[i], 0.25);
  }

}



}

