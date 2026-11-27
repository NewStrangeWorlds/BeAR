// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <iomanip>

#include "atmosphere.h"

#include "../../chemistry/chem_species.h"
#include "../../additional/physical_const.h"
#include "../../additional/exceptions.h"
#include "../../CUDA_kernels/data_management_kernels.h"
#include "../../chemistry/fixed_chemistry.h"


namespace bear{


Atmosphere::Atmosphere(
  const size_t nb_grid_points_,
  const std::vector<double>& atmos_boundaries,
  const bool use_gpu) : nb_grid_points(nb_grid_points_)
{
  createPressureGrid(atmos_boundaries);

  //initialise temperatures, altitudes, and number_densities to 0
  temperature.assign(nb_grid_points, 0.0);
  altitude.assign(nb_grid_points, 0.0);
  scale_height.assign(nb_grid_points, 0.0);
  number_densities.assign(
    nb_grid_points, std::vector<double>(constants::species_data.size(), 0.0));
  
  if (use_gpu)
  {
    allocateOnDevice(temperature_dev, nb_grid_points);
    allocateOnDevice(altitude_dev, nb_grid_points);
    allocateOnDevice(pressure_dev, nb_grid_points);
  }
}


bool Atmosphere::calcAtmosphereStructure(
  const double surface_gravity,
  const double reference_radius,
  const bool use_variable_gravity,
  Temperature* temperature_profile,
  const std::vector<double>& temp_parameters,
  std::vector<std::unique_ptr<Chemistry>>& chemistry,
  const std::vector<double>& chem_parameters_all)
{ 
  bool neglect_model = false;

  //temperature profile
  bool neglect_temperature = temperature_profile->calcProfile(
    temp_parameters, surface_gravity, pressure, temperature);

  if (neglect_temperature) neglect_model = true;


  //chemical composition
  std::vector<double> mean_molecular_weights(nb_grid_points, 0.0);

  number_densities.assign(
    nb_grid_points,
    std::vector<double>(constants::species_data.size(), 0.0));

  size_t nb_chem_param = 0;

  for (auto & i : chemistry)
  {
    std::vector<double> chem_parameters(
      chem_parameters_all.begin() + nb_chem_param,
      chem_parameters_all.begin() + nb_chem_param + i->nbParameters());
    
    nb_chem_param += i->nbParameters();
    
    bool neglect = i->calcChemicalComposition(
      chem_parameters, temperature, pressure, number_densities, mean_molecular_weights);
    
    if (neglect) neglect_model = true;
  }
  
  
  if (use_variable_gravity)
  {
    const bool unbound_atmosphere = calcAltitudeVariableGravity(
      surface_gravity, 
      reference_radius,
      mean_molecular_weights);

    if (unbound_atmosphere) neglect_model = true;
  }
  else
    calcAltitude(surface_gravity, reference_radius, mean_molecular_weights);
  
  calcScaleHeight(surface_gravity, reference_radius, use_variable_gravity, mean_molecular_weights);

  if (altitude_dev != nullptr)
  {
    moveToDevice(altitude_dev, altitude);
    moveToDevice(temperature_dev, temperature);
  }

  return neglect_model;
}



bool Atmosphere::calcAtmosphereStructure(
  const double surface_gravity,
  const double reference_radius,
  const bool use_variable_gravity,
  Temperature* temperature_profile,
  const std::vector<double>& temp_parameters,
  std::vector<std::unique_ptr<Chemistry>>& chemistry,
  const std::vector<double>& chem_parameters_all,
  const double mean_molecular_weight)
{ 
  bool neglect_model = false;

  //temperature profile
  bool neglect_temperature = temperature_profile->calcProfile(
    temp_parameters, surface_gravity, pressure, temperature);

  if (neglect_temperature) neglect_model = true;


  //chemical composition
  std::vector<double> mean_molecular_weights(nb_grid_points, 0.0);

  number_densities.assign(
    nb_grid_points,
    std::vector<double>(constants::species_data.size(), 0.0));

  size_t nb_chem_param = 0;

  for (auto & i : chemistry)
  {
    std::vector<double> chem_parameters(
      chem_parameters_all.begin() + nb_chem_param,
      chem_parameters_all.begin() + nb_chem_param + i->nbParameters());
    
    nb_chem_param += i->nbParameters();
    
    bool neglect = i->calcChemicalComposition(
      chem_parameters, temperature, pressure, number_densities, mean_molecular_weights);
    
    if (neglect) neglect_model = true;
  }
  
  mean_molecular_weights.assign(nb_grid_points, mean_molecular_weight);
  
  if (use_variable_gravity)
  {
    const bool unbound_atmosphere = calcAltitudeVariableGravity(
      surface_gravity, 
      reference_radius,
      mean_molecular_weights);

    if (unbound_atmosphere) neglect_model = true;
  }
  else
    calcAltitude(surface_gravity, reference_radius, mean_molecular_weights);
  
  calcScaleHeight(surface_gravity, reference_radius, use_variable_gravity, mean_molecular_weights);

  if (altitude_dev != nullptr)
  {
    moveToDevice(altitude_dev, altitude);
    moveToDevice(temperature_dev, temperature);
  }

  return neglect_model;
}



bool Atmosphere::calcAtmosphereStructure(
  const double surface_gravity,
  const double reference_radius,
  const double constant_scale_height,
  Temperature* temperature_profile,
  const std::vector<double>& temp_parameters,
  std::vector<std::unique_ptr<Chemistry>>& chemistry,
  const std::vector<double>& chem_parameters_all)
{ 
  bool neglect_model = false;

  //temperature profile
  bool neglect_temperature = temperature_profile->calcProfile(
    temp_parameters, surface_gravity, pressure, temperature);

  if (neglect_temperature) neglect_model = true;


  //chemical composition
  std::vector<double> mean_molecular_weights(nb_grid_points, 0.0);

  number_densities.assign(
    nb_grid_points,
    std::vector<double>(constants::species_data.size(), 0.0));

  size_t nb_chem_param = 0;

  for (auto & i : chemistry)
  {
    std::vector<double> chem_parameters(
      chem_parameters_all.begin() + nb_chem_param,
      chem_parameters_all.begin() + nb_chem_param + i->nbParameters());
    
    nb_chem_param += i->nbParameters();
    
    bool neglect = i->calcChemicalComposition(
      chem_parameters, temperature, pressure, number_densities, mean_molecular_weights);
    
    if (neglect) neglect_model = true;
  }

  calcAltitude(constant_scale_height, reference_radius);
  scale_height.assign(nb_grid_points, constant_scale_height*100000.0);

  if (altitude_dev != nullptr)
  {
    moveToDevice(altitude_dev, altitude);
    moveToDevice(temperature_dev, temperature);
  }

  return neglect_model;
}



void Atmosphere::setAtmosphericStructure(
  const double surface_gravity,
  const double reference_radius,
  const bool use_variable_gravity,
  const std::vector<double>& pressure_,
  const std::vector<double>& temperature_,
  const std::vector<std::string>& species_symbols,
  const std::vector< std::vector<double>>& mixing_ratios)
{ 
  temperature = temperature_;
  pressure = pressure_;

  FixedChemistry fixed_chemistry = FixedChemistry(species_symbols);
  
  std::vector<double> mean_molecular_weights(nb_grid_points, 0.0);

  number_densities.assign(
    nb_grid_points,
    std::vector<double>(constants::species_data.size(), 0.0));
  
  fixed_chemistry.setChemicalComposition(
    pressure, temperature, mixing_ratios, number_densities, mean_molecular_weights);

  if (use_variable_gravity)
    calcAltitudeVariableGravity(
      surface_gravity, 
      reference_radius,
      mean_molecular_weights);
  else
    calcAltitude(surface_gravity, reference_radius, mean_molecular_weights);
  
  calcScaleHeight(surface_gravity, reference_radius, use_variable_gravity, mean_molecular_weights);

  if (altitude_dev != nullptr)
  {
    moveToDevice(altitude_dev, altitude);
    moveToDevice(temperature_dev, temperature);
  }
}



//returns the pressure of the reference level, i.e. the level where
//the reference radius and the surface gravity are defined
//if no reference pressure has been set, the bottom of the grid is used
double Atmosphere::referencePressure() const
{
  if (reference_pressure <= 0.0) return pressure.front();

  //small tolerance since the grid boundaries are affected by round-off errors
  const double tolerance = 1e-10;

  if (reference_pressure > pressure.front()*(1.0 + tolerance) 
    || reference_pressure < pressure.back()*(1.0 - tolerance))
  {
    std::ostringstream message;
    message << "The reference pressure of " << reference_pressure 
            << " bar is outside of the pressure grid ("
            << pressure.front() << " to " << pressure.back() << " bar).\n";

    throw InvalidInput(std::string ("Atmosphere::referencePressure"), message.str());
  }

  return std::min(std::max(reference_pressure, pressure.back()), pressure.front());
}


//integrates a level quantity f over ln(p) with the trapezoidal rule
//returns the integral of f dln(p) between the reference pressure and each level
//it is positive above the reference pressure and negative below it
//if the reference pressure is located between two levels, f is 
//interpolated linearly in ln(p)
std::vector<double> Atmosphere::integrateFromReferencePressure(
  const std::vector<double>& f) const
{
  std::vector<double> integral(nb_grid_points, 0.0);

  if (nb_grid_points < 2) return integral;

  const double p_ref = referencePressure();

  //index of the level directly below the reference pressure
  size_t k = 0;

  while (k < nb_grid_points-2 && pressure[k+1] >= p_ref) ++k;

  const double weight = std::log(pressure[k]/p_ref) / std::log(pressure[k]/pressure[k+1]);
  const double f_ref = f[k] + weight * (f[k+1] - f[k]);

  //the two levels enclosing the reference pressure
  integral[k+1] = 0.5 * (f_ref + f[k+1]) * std::log(p_ref/pressure[k+1]);
  integral[k] = - 0.5 * (f_ref + f[k]) * std::log(pressure[k]/p_ref);

  //upwards
  for (size_t i=k+2; i<nb_grid_points; ++i)
    integral[i] = integral[i-1] 
                + 0.5 * (f[i] + f[i-1]) * std::log(pressure[i-1]/pressure[i]);

  //downwards
  for (size_t i=k; i-- > 0; )
    integral[i] = integral[i+1] 
                - 0.5 * (f[i] + f[i+1]) * std::log(pressure[i]/pressure[i+1]);

  return integral;
}


//determine the vertical grid via hydrostatic equilibrium
//integrates dz = -H dln(p) with the trapezoidal rule in ln(p)
//this is exact for a layer with a constant scale height
//the reference radius is the radius at the reference pressure
//the altitudes are given with respect to the bottom of the grid
void Atmosphere::calcAltitude(
  const double surface_gravity, 
  const double reference_radius,
  const std::vector<double>& mean_molecular_weights)
{
  std::vector<double> pressure_scale_height(nb_grid_points, 0.0);

  for (size_t i=0; i<nb_grid_points; ++i)
    pressure_scale_height[i] = constants::gas_constant * temperature[i] 
                             / (mean_molecular_weights[i] * surface_gravity);

  //altitudes with respect to the reference pressure
  altitude = integrateFromReferencePressure(pressure_scale_height);

  bottom_radius = reference_radius + altitude.front();

  const double bottom_altitude = altitude.front();

  for (auto & z : altitude) z -= bottom_altitude;
}


//determine the vertical grid via hydrostatic equilibrium
//with a gravity that varies as g(r) = g_ref * (r_ref/r)^2
//where g_ref and r_ref are the values at the reference pressure
//in this case, hydrostatic equilibrium can be written as
//d(1/r) = R*T/(mu * g_ref * r_ref^2) dln(p)
//which is integrated with the trapezoidal rule in ln(p)
//this is exact for a layer with constant T/mu
//the altitudes are given with respect to the bottom of the grid
//returns true if the atmosphere is not gravitationally bound
bool Atmosphere::calcAltitudeVariableGravity(
  const double surface_gravity_reference,
  const double reference_radius,
  const std::vector<double>& mean_molecular_weights)
{ 
  altitude.assign(nb_grid_points, 0.0);
  bottom_radius = reference_radius;

  std::vector<double> integrand(nb_grid_points, 0.0);

  for (size_t i=0; i<nb_grid_points; ++i)
    integrand[i] = constants::gas_constant * temperature[i] 
                 / (mean_molecular_weights[i] * surface_gravity_reference * reference_radius*reference_radius);

  //1/r_ref - 1/r
  const std::vector<double> delta_inverse_radius = integrateFromReferencePressure(integrand);

  std::vector<double> radius(nb_grid_points, 0.0);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    const double inverse_radius = 1.0/reference_radius - delta_inverse_radius[i];

    //the pressure never drops to the value of this level at any finite radius
    if (inverse_radius <= 0.0) return true;

    radius[i] = 1.0/inverse_radius;
  }

  bottom_radius = radius.front();

  for (size_t i=0; i<nb_grid_points; ++i)
    altitude[i] = radius[i] - bottom_radius;

  return false;
}


//determine the vertical grid for a constant scale height
void Atmosphere::calcAltitude(
  const double constant_scale_height,
  const double reference_radius)
{
  altitude.assign(nb_grid_points, 0.0);
  
  //from km to cm
  const double scale_height = constant_scale_height * 100000.0;

  for (size_t i=1; i<nb_grid_points; i++)
  {
    altitude[i] = - scale_height * std::log(pressure[i]/pressure[0]);
  }

  bottom_radius = reference_radius - scale_height * std::log(pressure[0]/referencePressure());
}


//pressure scale height at each level
//for a variable gravity, the local gravity g(r) = g_ref * (r_ref/r)^2 is used
void Atmosphere::calcScaleHeight(
  const double surface_gravity, 
  const double reference_radius,
  const bool use_variable_gravity,
  const std::vector<double>& mean_molecular_weights)
{
  scale_height.assign(nb_grid_points, 0.0);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    double gravity = surface_gravity;

    if (use_variable_gravity)
    {
      const double radius = bottom_radius + altitude[i];
      gravity *= (reference_radius/radius) * (reference_radius/radius);
    }

    scale_height[i] = constants::gas_constant * temperature[i] / (mean_molecular_weights[i] * gravity);
  }
}


//Creates a pressure grid between the two boundary points
//nb_grid_point pressures vales are equidistantly spread in the log(p) space
void Atmosphere::createPressureGrid(const std::vector<double>& atmos_boundaries)
{
  pressure.assign(nb_grid_points, 0.0);

  const double min_pressure = atmos_boundaries[1];
  const double max_pressure = atmos_boundaries[0];

  pressure.front() = std::log10(max_pressure);

  const double log_step = (std::log10(max_pressure) - std::log10(min_pressure)) / (nb_grid_points - 1.0);

  for (size_t i=1; i<nb_grid_points-1; ++i)
    pressure[i] = pressure[i-1] - log_step;

  for (size_t i=0; i<nb_grid_points-1; ++i)
    pressure[i] = pow(10.0, pressure[i]);

  pressure.back() = min_pressure;
}


Atmosphere::~Atmosphere()
{
  if (temperature_dev != nullptr) deleteFromDevice(temperature_dev);
  if (altitude_dev != nullptr) deleteFromDevice(altitude_dev);
  if (pressure_dev != nullptr) deleteFromDevice(pressure_dev);
}


}
