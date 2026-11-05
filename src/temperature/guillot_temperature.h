// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _guillot_temperature_h
#define _guillot_temperature_h

#include <iostream>
#include <vector>
#include <string>

#include "temperature.h"


namespace bear {


class GuillotTemperature : public Temperature{
  public:
    GuillotTemperature(const std::string profile_type);
    virtual ~GuillotTemperature() {}
    virtual bool calcProfile(
      const std::vector<double>& parameters,
      const double surface_gravity,
      const std::vector<double>& pressure,
      std::vector<double>& temperature);
  private:
    unsigned int profile = 0;

    void profileBeamSource(
      const std::vector<double>& optical_depth,
      const double temperature_irr,
      const double temperature_int,
      const double mu,
      const double gamma,
      std::vector<double>& temperature);
    void profileIsotropicSource(
      const std::vector<double>& optical_depth,
      const double temperature_irr,
      const double temperature_int,
      const double flux_distribution,
      const double gamma,
      std::vector<double>& temperature);
};


}
#endif 
