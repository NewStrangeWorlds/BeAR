// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _cubic_b_spline_temperature_h
#define _cubic_b_spline_temperature_h

#include "temperature.h"

#include <vector>


namespace bear {


class CubicBSplineTemperature : public Temperature{
  public:
    CubicBSplineTemperature(
      const size_t nb_control_points_,
      const Parametrisation parametrisation_);
    virtual ~CubicBSplineTemperature() {}
    virtual bool calcProfile(
      const std::vector<double>& parameters,
      const double surface_gravity,
      const std::vector<double>& pressure,
      std::vector<double>& temperature);
  private:
     const size_t nb_control_points;
};


}
#endif 
