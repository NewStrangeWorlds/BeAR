
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


#include <iostream>
#include <omp.h>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "transport_coeff.h"

#include "species_definition.h"
#include "../spectral_grid/spectral_grid.h"
#include "../chemistry/chem_species.h"
#include "../config/global_config.h"
#include "../additional/exceptions.h"


namespace bear{


TransportCoefficients::TransportCoefficients(
  GlobalConfig* config_ptr,
  SpectralGrid* grid_ptr, 
  const std::vector<std::string>& opacity_species_symbol,
  const std::vector<std::string>& opacity_species_folder)
{
  config = config_ptr;
  spectral_grid = grid_ptr;


  std::vector<size_t> spectral_indices = spectral_grid->spectralIndexList();
  
  gas_species.reserve(opacity_species_symbol.size());

  bool all_species_added = true;

  for (size_t i=0; i<opacity_species_symbol.size(); ++i)
  {
    bool added = addOpacitySpecies(opacity_species_symbol[i], opacity_species_folder[i]);

    if (!added) all_species_added = false;
  }

  
  std::cout << "\nOpacity specied added:\n";
  for (auto & i : gas_species)
  {
    std::string data_available = "data available: No";
    
    if (i->dataAvailable())
      data_available = "data available: Yes";

    std::cout  << std::setw(8) << std::left << i->species_name << "\t" 
               << std::setw(40) << std::left << i->species_folder << "\t" 
               << std::setw(10) << std::left << data_available << "\n\n";
  }
  
  
  for (auto & i : gas_species)
  {
    if (i->dataAvailable() == false)
    {
      std::string error_message = "Unable to locate all opacity data.\n";
      throw InvalidInput(std::string ("TransportCoefficients::TransportCoefficients"), error_message);
    }
  }

  if (!all_species_added)
  {
    std::string error_message = "Not all opacities species from the model config file could be added!\nThis could be caused by a species unknown to BeAR.\n";
    throw InvalidInput(std::string ("TransportCoefficients::TransportCoefficients"), error_message);
  }

}




bool TransportCoefficients::addOpacitySpecies(
  const std::string& species_symbol, const std::string& species_folder)
{
  //first, the species cases for which separate classes are available
  if (species_symbol == "CIA-H2-H2")
  {
    gas_species.push_back(
      std::make_unique<GasGeneric>(
        config,
        spectral_grid,
        _H2,
        "CIA H2-H2",
        species_folder,
        std::vector<size_t>{_H2}));

    return true;
  }

  if (species_symbol == "CIA-H2-He")
  {
    gas_species.push_back(
      std::make_unique<GasGeneric>(
        config,
        spectral_grid,
        _H2,
        "CIA H2-He",
        species_folder,
        std::vector<size_t>{_He}));

    return true;
  }


  if (species_symbol == "CIA-H-He")
  {
    gas_species.push_back(
      std::make_unique<GasGeneric>(
        config,
        spectral_grid,
        _H,
        "CIA H-He",
        species_folder,
        std::vector<size_t>{_He}));

    return true;
  }

  //H- free-free and bound-free continuum
  if (species_symbol == "H-")
  {
    gas_species.push_back(std::make_unique<GasHm>(config, spectral_grid));

    return true;
  }

  //H2 Rayleigh scattering
  if (species_symbol == "H2" && species_folder == "Rayleigh")
  {
    gas_species.push_back(std::make_unique<GasH2Rayleigh>(config, spectral_grid, ""));

    return true;
  }

  //He Rayleigh scattering
  if (species_symbol == "He" && species_folder == "Rayleigh")
  {
    gas_species.push_back(std::make_unique<GasHeRayleigh>(config, spectral_grid, ""));

    return true;
  }

  //H Rayleigh scattering
  if (species_symbol == "H" && species_folder == "Rayleigh")
  {
    gas_species.push_back(std::make_unique<GasHRayleigh>(config, spectral_grid, ""));

    return true;
  }

  //CO Rayleigh
  if (species_symbol == "CO" && species_folder == "Rayleigh")
  {
    gas_species.push_back(std::make_unique<GasCORayleigh>(config, spectral_grid, ""));

    return true;
  }


  //CO2 Rayleigh
  if (species_symbol == "CO2" && species_folder == "Rayleigh")
  {
    gas_species.push_back(std::make_unique<GasCO2Rayleigh>(config, spectral_grid, ""));

    return true;
  }


  //CH4 Rayleigh
  if (species_symbol == "CH4" && species_folder == "Rayleigh")
  {
    gas_species.push_back(std::make_unique<GasCH4Rayleigh>(config, spectral_grid, ""));

    return true;
  }


  //H2O Rayleigh
  if (species_symbol == "H2O" && species_folder == "Rayleigh")
  {
    gas_species.push_back(std::make_unique<GasH2ORayleigh>(config, spectral_grid, ""));

    return true;
  }


  //now we try the generic ones
  for (size_t i=0; i<constants::species_data.size(); ++i)
  {
    if (constants::species_data[i].symbol == species_symbol)
    {
      gas_species.push_back(
        std::make_unique<GasGeneric>(
          config,
          spectral_grid,
          constants::species_data[i].id,
          constants::species_data[i].symbol,
          species_folder));

      return true;
    }
  }

  //we haven't found the corresponding species
  std::cout << "Opacity species " 
    << species_symbol 
    << " has not been found in the internal list located in chem_species.h!\n";


  return false;
}


//calculates the transport coefficients on the CPU
//calls the calculation method of the individual opacity species
void TransportCoefficients::calculate(
  const double temperature,
  const double pressure,
  const std::vector<double>& number_densities,
  std::vector<double>& absorption_coeff,
  std::vector<double>& scattering_coeff)
{
  absorption_coeff.assign(spectral_grid->nbSpectralPoints(), 0);
  scattering_coeff.assign(spectral_grid->nbSpectralPoints(), 0);


  for (unsigned int i=0; i<gas_species.size(); i++)
    gas_species[i]->calcTransportCoefficients(
      temperature,
      pressure,
      number_densities,
      absorption_coeff,
      scattering_coeff);
}



//calculates the transport coefficients on the GPU
//calculations are stored on the GPU, nothing is returned
//the layer coefficients are a temporary storage for a given p-T point
void TransportCoefficients::calculateGPU(
  const double temperature,
  const double pressure,
  const std::vector<double>& number_densities,
  const size_t nb_grid_points,
  const size_t grid_point,
  float* absorption_coeff_device,
  float* scattering_coeff_device)
{
  for (unsigned int i=0; i<gas_species.size(); i++)
    gas_species[i]->calcTransportCoefficientsGPU(
      temperature,
      pressure,
      number_densities,
      nb_grid_points,
      grid_point,
      absorption_coeff_device,
      scattering_coeff_device);
}



void TransportCoefficients::prepareBatchedGPU(
  const Atmosphere& atmosphere,
  std::vector<float*>& cs1_ptrs, std::vector<float*>& cs2_ptrs,
  std::vector<float*>& cs3_ptrs, std::vector<float*>& cs4_ptrs,
  std::vector<float>& temp_factors, std::vector<float>& pres_factors,
  std::vector<float>& cs_log_number_densities, std::vector<int>& cs_grid_points,
  std::vector<float*>& ray_ptrs,
  std::vector<double>& ray_number_densities, std::vector<int>& ray_grid_points)
{
  cs1_ptrs.clear(); cs2_ptrs.clear(); cs3_ptrs.clear(); cs4_ptrs.clear();
  temp_factors.clear(); pres_factors.clear();
  cs_log_number_densities.clear(); cs_grid_points.clear();
  ray_ptrs.clear(); ray_number_densities.clear(); ray_grid_points.clear();

  const size_t nb_grid_points = atmosphere.nb_grid_points;

  for (size_t gp = 0; gp < nb_grid_points; ++gp)
  {
    const double temperature = atmosphere.temperature[gp];
    const double pressure = atmosphere.pressure[gp];
    const auto& number_densities = atmosphere.number_densities[gp];

    for (auto& species : gas_species)
    {
      double number_density = number_densities[species->species_index];

      for (const auto& partner : species->getCIACollisionPartners())
        number_density *= number_densities[partner];

      if (number_density == 0) continue;

      double reference_pressure = pressure;
      if (species->getPressureReferenceSpecies() != _TOTAL)
        reference_pressure *= number_densities[species->getPressureReferenceSpecies()]
                            / number_densities[_TOTAL];

      if (species->hasCrossSections())
      {
        auto meta = species->prepareCrossSectionMetadata(reference_pressure, temperature);

        if (meta.valid)
        {
          cs1_ptrs.push_back(meta.cs1);
          cs2_ptrs.push_back(meta.cs2);
          cs3_ptrs.push_back(meta.cs3);
          cs4_ptrs.push_back(meta.cs4);
          temp_factors.push_back(meta.temperature_interpol_factor);
          pres_factors.push_back(meta.pressure_interpol_factor);
          cs_log_number_densities.push_back(static_cast<float>(std::log10(number_density)));
          cs_grid_points.push_back(static_cast<int>(gp));
        }
      }

      if (species->hasRayleigh())
      {
        float* ray_dev = species->getRayleighDevicePtr();
        if (ray_dev != nullptr)
        {
          ray_ptrs.push_back(ray_dev);
          ray_number_densities.push_back(number_densities[species->species_index]);
          ray_grid_points.push_back(static_cast<int>(gp));
        }
      }
    }
  }
}


void TransportCoefficients::calculateContinuumGPU(
  const double temperature,
  const double pressure,
  const std::vector<double>& number_densities,
  const size_t nb_grid_points,
  const size_t grid_point,
  float* absorption_coeff_device)
{
  for (auto& species : gas_species)
    species->calcContinuumGPU(
      temperature, number_densities, nb_grid_points, grid_point, absorption_coeff_device);
}


TransportCoefficients::~TransportCoefficients()
{
}



}
