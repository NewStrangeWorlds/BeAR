// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

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
      const std::vector<double>& atmos_boundaries,
      const Parametrisation parametrisation_);
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
