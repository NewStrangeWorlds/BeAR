/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2020 Daniel Kitzmann
*
* Helios-r2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Helios-r2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* Helios-r2 directory under <LICENSE>. If not, see
* <http://www.gnu.org/licenses/>.
*/


#ifndef _piecewise_poly_temperature_h
#define _piecewise_poly_temperature_h

#include "temperature.h"

#include "../additional/piecewise_poly.h"

#include <vector>



namespace helios {


class PiecewisePolynomialTemperature : public Temperature{
  public:
    PiecewisePolynomialTemperature(const size_t nb_elements_in, const size_t polynomial_degree_in, const double atmos_boundaries [2]);
    virtual ~PiecewisePolynomialTemperature() {}
    virtual bool calcProfile(const std::vector<double>& parameters, const std::vector<double>& pressure, std::vector<double>& temperature);
  private:
    PiecewisePolynomial temperature_profile;
    const size_t nb_elements {}; 
    const size_t polynomial_degree {};
};


}
#endif 
