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


#ifndef _madhusudhan_seager_temperature_h
#define _madhusudhan_seager_temperature_h

#include <vector>
#include <string>

#include "temperature.h"


namespace bear {


// Madhusudhan & Seager (2009) three-layer T-P profile.
//
// Free parameters (6):
//   T0      - temperature at the top of the atmosphere [K]
//   logP1   - log10(P1 / bar): boundary between layers 1 and 2
//   logP2   - log10(P2 / bar): reference pressure for layer 2 (temperature minimum)
//   logP3   - log10(P3 / bar): pressure below which the atmosphere is isothermal
//   alpha1  - gradient parameter for layer 1  (> 0)
//   alpha2  - gradient parameter for layer 2  (> 0)
//
// Profile (P0 = pressure at top of atmosphere = pressure.back()):
//   Layer 1 (P < P1):           T(P) = [ln(P/P0)/alpha1]^2 + T0
//   Layer 2 (P1 <= P < P3):     T(P) = [ln(P/P2)/alpha2]^2 + T2
//     where T2 = [ln(P1/P0)/alpha1]^2 + T0 - [ln(P1/P2)/alpha2]^2  (continuity at P1)
//   Layer 3 (P >= P3):          T = T(P3)  [isothermal]
class MadhusudhanSeagerTemperature : public Temperature {
  public:
    MadhusudhanSeagerTemperature();
    virtual ~MadhusudhanSeagerTemperature() {}
    virtual bool calcProfile(
      const std::vector<double>& parameters,
      const double surface_gravity,
      const std::vector<double>& pressure,
      std::vector<double>& temperature);
};


}
#endif
