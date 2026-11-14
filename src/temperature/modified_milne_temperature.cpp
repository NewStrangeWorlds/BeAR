// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include "modified_milne_temperature.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"

#include <algorithm>
#include <vector>
#include <cmath>


namespace bear {


//calculate the temperature based on the modified Milne solution
//  T^4(p) = 3/4 T_eff^4 (q(p) + tau(p))
//  q(p)   = q_inf + (q_0 - q_inf) exp[-(p/p_t)^beta]
//  q_0    = 4/3 (T_skin/T_eff)^4
bool ModifiedMilneTemperature::calcProfile(
  const std::vector<double>& parameters, 
  const double surface_gravity,
  const std::vector<double>& pressure, 
  std::vector<double>& temperature)
{
  temperature.assign(pressure.size(), 0);

  const double kappa_ross = parameters[0];
  const double t_eff = parameters[1];
  const double t_ratio = parameters[2];         //T_skin / T_eff
  const double p_transition = parameters[3];    //in bar
  const double beta = parameters[4];

  const double q_0 = 4.0/3.0 * std::pow(t_ratio, 4);

  for (size_t i=0; i<temperature.size(); ++i)
  {
    const double tau_ross = kappa_ross * pressure[i] * 1e6 / surface_gravity;

    const double sigmoid = std::exp(-std::pow(pressure[i] / p_transition, beta));
    const double q = q_inf + (q_0 - q_inf) * sigmoid;

    temperature[i] = 3.0/4.0 * std::pow(t_eff, 4) * (q + tau_ross);
    temperature[i] = std::pow(temperature[i], 0.25);
  }

  //neglect models with too low temperatures
  bool neglect_model = false;
  
  for (auto & i : temperature)
    if (i < 50) {i = 50; neglect_model = true;}


  return neglect_model;
}


}
