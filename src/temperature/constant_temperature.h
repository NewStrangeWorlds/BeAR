// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _constant_temperature_h
#define _constant_temperature_h

#include "temperature.h"


#include <vector>



namespace bear {


class ConstantTemperature : public Temperature{
  public:
    ConstantTemperature() {parameter_names = {"temperature"};}
    virtual ~ConstantTemperature() {}
    virtual bool calcProfile(
      const std::vector<double>& parameters,
      const double surface_gravity,
      const std::vector<double>& pressure,
      std::vector<double>& temperature) {
        temperature.assign(pressure.size(), parameters[0]);
        return false;
      }
};


}
#endif 
