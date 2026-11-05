// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _milne_solution_temperature_h
#define _milne_solution_temperature_h

#include <iostream>
#include <vector>

#include "temperature.h"


namespace bear {


class MilneTemperature : public Temperature{
  public:
    MilneTemperature() {
        parameter_names = {"temp_kappa_ross", "temp_eff"};
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
