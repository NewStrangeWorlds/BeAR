// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _modified_milne_temperature_h
#define _modified_milne_temperature_h

#include <iostream>
#include <vector>

#include "temperature.h"


namespace bear {


//Milne's solution with a modified Hopf function that transitions via a
//stretched exponential in pressure from a free skin temperature at zero
//optical depth to the classic Milne solution at depth
//(follows the "Modified Milne" profile of the Exo_Skryer code by E. Lee)
class ModifiedMilneTemperature : public Temperature{
  public:
    ModifiedMilneTemperature() {
        parameter_names = {
          "temp_kappa_ross", "temp_eff", "temp_ratio", "temp_p_t", "temp_beta"};
        std::cout << "\n- Temperature profile: modified Milne's solution\n\n";
      }
    virtual ~ModifiedMilneTemperature() {}
    virtual bool calcProfile(
      const std::vector<double>& parameters,
      const double surface_gravity,
      const std::vector<double>& pressure,
      std::vector<double>& temperature);
  private:
    //asymptotic value of the Hopf function for tau -> infinity
    const double q_inf = 0.710446;
};


}
#endif 
