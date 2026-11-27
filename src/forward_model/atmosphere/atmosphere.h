// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _atmosphere_h
#define _atmosphere_h

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <memory>

#include "../../chemistry/chemistry.h"
#include "../../temperature/temperature.h"


namespace bear {


class Atmosphere {
  public:
    Atmosphere (
      const size_t nb_grid_points_,
      const std::vector<double>& atmos_boundaries,
      const bool use_gpu);
    ~Atmosphere();

    const size_t nb_grid_points = 0;

    std::vector<double> pressure;
    std::vector<double> temperature;
    std::vector<double> altitude;
    std::vector<double> scale_height;
    std::vector< std::vector<double> > number_densities;

    //pressure in bar where the reference radius and the surface gravity are defined
    //for values <= 0, the bottom of the grid is used
    double reference_pressure = 0.0;

    //radius of the lowest grid point, the altitudes are given with respect to it
    //calculated together with the altitudes from the reference radius
    double bottom_radius = 0.0;

    float* altitude_dev = nullptr;
    float* pressure_dev = nullptr;
    float* temperature_dev = nullptr;

    bool calcAtmosphereStructure(
      const double surface_gravity,
      const double reference_radius,
      const bool use_variable_gravity,
      Temperature* temperature_profile,
      const std::vector<double>& temp_parameters,
      std::vector<std::unique_ptr<Chemistry>>& chemistry,
      const std::vector<double>& chem_parameters);
    bool calcAtmosphereStructure(
      const double surface_gravity,
      const double reference_radius,
      const bool use_variable_gravity,
      Temperature* temperature_profile,
      const std::vector<double>& temp_parameters,
      std::vector<std::unique_ptr<Chemistry>>& chemistry,
      const std::vector<double>& chem_parameters,
      const double mean_molecular_weight);
    bool calcAtmosphereStructure(
      const double surface_gravity,
      const double reference_radius,
      const double constant_scale_height,
      Temperature* temperature_profile,
      const std::vector<double>& temp_parameters,
      std::vector<std::unique_ptr<Chemistry>>& chemistry,
      const std::vector<double>& chem_parameters);

    void setAtmosphericStructure(
      const double surface_gravity,
      const double reference_radius,
      const bool use_variable_gravity,
      const std::vector<double>& pressure_,
      const std::vector<double>& temperature_,
      const std::vector<std::string>& species_symbols,
      const std::vector< std::vector<double>>& mixing_ratios);

  private:
    void createPressureGrid(const std::vector<double>& domain_boundaries);
    double referencePressure() const;
    std::vector<double> integrateFromReferencePressure(
      const std::vector<double>& f) const;
    void calcAltitude(
      const double g, 
      const double reference_radius,
      const std::vector<double>& mean_molecular_weights);
    bool calcAltitudeVariableGravity(
      const double g, 
      const double reference_radius,
      const std::vector<double>& mean_molecular_weights);
    void calcAltitude(
      const double constant_scale_height,
      const double reference_radius);
    void calcScaleHeight(
      const double g, const std::vector<double>& mean_molecular_weights);
};


}


#endif
