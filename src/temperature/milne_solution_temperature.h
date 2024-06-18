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


#ifndef _milne_solution_temperature_h
#define _milne_solution_temperature_h

#include <iostream>
#include <vector>

#include "temperature.h"


namespace helios {


class MilneTemperature : public Temperature{
  public:
    MilneTemperature() {
        nb_parameters = 2;
        std::cout << "\n- Temperature profile: Milne's solution\n\n";
      }
    virtual ~MilneTemperature() {}
    virtual bool calcProfile(
      const std::vector<double>& parameters,
      const double surface_gravity,
      const std::vector<double>& pressure,
      std::vector<double>& temperature);
  private:
    //fit coefficients for the Hopf function
    const std::vector<double> fit_p {0.6162, -0.3799, 2.395, -2.041, 2.578};
    const std::vector<double> fit_q {-0.9799, 3.917, -3.17, 3.69};

    double hopfFunction(const double optical_depth);
};


}
#endif 
