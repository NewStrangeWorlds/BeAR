// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

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
      std::vector<double>& temperature,
      const double pressure_rcb,
      const double temperature_rcb);

     const size_t nb_control_points;
};


}
#endif 
