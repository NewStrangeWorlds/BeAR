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


#ifndef _piecewise_poly_temperature_h
#define _piecewise_poly_temperature_h

#include "temperature.h"

#include "../additional/piecewise_poly.h"

#include <vector>



namespace bear {


class PiecewisePolynomialTemperature : public Temperature{
  public:
    PiecewisePolynomialTemperature(
      const size_t nb_elements_in,
      const size_t polynomial_degree_in,
      const std::vector<double>& atmos_boundaries);
    virtual ~PiecewisePolynomialTemperature() {}
    virtual bool calcProfile(
      const std::vector<double>& parameters,
      const double surface_gravity,
      const std::vector<double>& pressure,
      std::vector<double>& temperature);
  private:
    PiecewisePolynomial temperature_profile;
    const size_t nb_elements {}; 
    const size_t polynomial_degree {};
};


}
#endif 
