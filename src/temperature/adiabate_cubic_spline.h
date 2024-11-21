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


#ifndef _adiabate_spline_temperature_h
#define _adiabate_spline_temperature_h

#include "temperature.h"

#include <vector>


namespace bear {


class AdiabateSplineTemperature : public Temperature{
  public:
    AdiabateSplineTemperature(const size_t nb_control_points_);
    virtual ~AdiabateSplineTemperature() {}
    virtual bool calcProfile(
      const std::vector<double>& parameters,
      const double surface_gravity,
      const std::vector<double>& pressure,
      std::vector<double>& temperature);
  private:
     int findRadiavativeConvectiveBoundary(
      const double rcb_pressure,
      const std::vector<double>& pressure);
     void addAdiabate(
      const double temperature_bottom,
      const unsigned int rcb_idx,
      const double gamma,
      const std::vector<double>& pressure,
      std::vector<double>& temperature);
    void addCubicSpline(
      const std::vector<double>& parameters,
      const unsigned int rcb_idx,
      const std::vector<double>& pressure,
      std::vector<double>& temperature);

     const size_t nb_control_points;
};


}
#endif 
