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


#ifndef _atmosphere_h
#define _atmosphere_h

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

#include "../../chemistry/chemistry.h"
#include "../../temperature/temperature.h"


namespace helios {


class Atmosphere {
  public:
    Atmosphere (
      const size_t nb_grid_points_,
      const double atmos_boundaries [2],
      const bool use_gpu);
    ~Atmosphere();

    const size_t nb_grid_points = 0;

    std::vector<double> pressure;
    std::vector<double> temperature;
    std::vector<double> altitude;
    std::vector<double> scale_height;
    std::vector< std::vector<double> > number_densities;

    double* altitude_dev = nullptr;
    double* pressure_dev = nullptr;
    double* temperature_dev = nullptr;

    bool calcAtmosphereStructure(
      const double surface_gravity,
      Temperature* temperature_profile,
      const std::vector<double>& temp_parameters,
      std::vector<Chemistry*>& chemistry,
      const std::vector<double>& chem_parameters);
  private:
    void createPressureGrid(const double domain_boundaries [2]);
    void calcAltitude(
      const double g, const std::vector<double>& mean_molecular_weights);
    void calcScaleHeight(
      const double g, const std::vector<double>& mean_molecular_weights);
};


}


#endif
