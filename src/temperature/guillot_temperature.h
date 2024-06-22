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
