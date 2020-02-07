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


#include "transport_coeff.h"


#include "species_definition.h"
#include "../chemistry/chem_species.h"
#include "transport_coeff_single_species.h"
#include "../spectral_grid/spectral_grid.h"
#include "../CUDA_kernels/data_management_kernels.h"


#include <iostream>
#include <omp.h>
#include <fstream>
#include <iomanip>
#include <cmath>


#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>



namespace helios{


TransportCoefficients::TransportCoefficients(GlobalConfig* config_ptr, SpectralGrid* grid_ptr, 
                                             const std::vector<std::string>& opacity_species_symbol, const std::vector<std::string>& opacity_species_folder)
{
  config = config_ptr;
  spectral_grid = grid_ptr;


  std::vector<size_t> spectral_indices = spectral_grid->spectralIndexList();


  bool all_species_added = true;

  for (size_t i=0; i<opacity_species_symbol.size(); ++i)
  {
    bool added = addOpacitySpecies(opacity_species_symbol[i], opacity_species_folder[i]);

    if (!added) all_species_added = false;
  }

  
  std::cout << "\nOpacity specied added:\n";
  for (auto & i : gas_species)
    std::cout << i->getSpeciesName() << "\n";
    
  if (!all_species_added) 
    std::cout << "\nWarning, not all opacities species from the model config file could be added!\n";


  for (auto & i : gas_species)
    i->init(spectral_indices);
}



bool TransportCoefficients::addOpacitySpecies(const std::string& species_symbol, const std::string& species_folder)
{
  //first, the species cases for which separate classes are available
  if (species_symbol == "H2")
  {
    GasH2* h2 = new GasH2(config, spectral_grid);
    gas_species.push_back(h2);
   
    return true;
  }


  //now we try the generic ones
  for (size_t i=0; i<constants::species_data.size(); ++i)
  {
    if (constants::species_data[i].symbol == species_symbol)
    {
      gas_species.push_back(new GasGeneric(config, spectral_grid, i, constants::species_data[i].symbol));
      
      return true;
    } 
  }
  
  
  //we haven't found the corresponding species
  std::cout << "Opacity species " << species_symbol << " has not been found in the internal list located in chem_species.h!\n";


  return false;
}




TransportCoefficients::~TransportCoefficients()
{

  for (unsigned int i=0; i<gas_species.size(); ++i)
    delete gas_species[i];

}



//calculates the transport coefficients on the CPU
//calls the calculation method of the individual opacity species
void TransportCoefficients::calcTransportCoefficients(const double temperature, const double pressure, const std::vector<double>& number_densities,
                                                      std::vector<double>& absorption_coeff, std::vector<double>& scattering_coeff)
{
  omp_set_nested(1);

  #pragma omp parallel for schedule(dynamic, 1)
  for (unsigned int i=0; i<gas_species.size(); i++)
    gas_species[i]->prepareCalculation(temperature, pressure);

  omp_set_nested(0);


  absorption_coeff.assign(spectral_grid->nbSpectralPoints(), 0);
  scattering_coeff.assign(spectral_grid->nbSpectralPoints(), 0);


  for (unsigned int i=0; i<gas_species.size(); i++)
    gas_species[i]->calcTransportCoefficients(temperature, pressure, number_densities, absorption_coeff, scattering_coeff);
}




//calculates the transport coefficients on the GPU
//calculations are stored on the GPU, nothing is returned
//the layer coefficients are a temporary storage for a given p-T point
void TransportCoefficients::calcTransportCoefficientsGPU(const double temperature, const double pressure, const std::vector<double>& number_densities,
                                                         const size_t nb_grid_points, const size_t grid_point,
                                                         double* absorption_coeff_device, double* scattering_coeff_device)
{
  omp_set_nested(1);

  #pragma omp parallel for schedule(dynamic, 1)
  for (unsigned int i=0; i<gas_species.size(); i++)
    gas_species[i]->prepareCalculation(temperature, pressure);

  omp_set_nested(0);


  for (unsigned int i=0; i<gas_species.size(); i++)
    gas_species[i]->calcTransportCoefficientsGPU(temperature, pressure, number_densities,
                                                 nb_grid_points, grid_point,
                                                 absorption_coeff_device, scattering_coeff_device);
}




}



