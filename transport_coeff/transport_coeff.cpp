/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2019 Daniel Kitzmann
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




void TransportCoefficients::init(GlobalConfig* config_ptr, SpectralGrid* grid_ptr)
{
  config = config_ptr;
  spectral_grid = grid_ptr;


  GasH2* h2 = new GasH2(config_ptr, spectral_grid);
  gas_species.push_back(h2);


  GasGeneric* na = new GasGeneric(config_ptr, spectral_grid, _Na, "Na");
  gas_species.push_back(na);


  GasGeneric* k = new GasGeneric(config_ptr, spectral_grid, _K, "K");
  gas_species.push_back(k);


  GasGeneric* h2o = new GasGeneric(config_ptr, spectral_grid, _H2O, "H2O");
  gas_species.push_back(h2o);


  GasGeneric* co = new GasGeneric(config_ptr, spectral_grid, _CO, "CO");
  gas_species.push_back(co);


  GasGeneric* ch4 = new GasGeneric(config_ptr, spectral_grid, _CH4, "CH4");
  gas_species.push_back(ch4);


  GasGeneric* co2 = new GasGeneric(config_ptr, spectral_grid, _CO2, "CO2");
  gas_species.push_back(co2);


  GasGeneric* nh3 = new GasGeneric(config_ptr, spectral_grid, _NH3, "NH3");
  gas_species.push_back(nh3);


  GasGeneric* h2s = new GasGeneric(config_ptr, spectral_grid, _H2S, "H2S");
  gas_species.push_back(h2s);


  std::vector<size_t> spectral_indices = spectral_grid->spectralIndexList();


  for (unsigned int i=0; i<gas_species.size(); ++i)
    std::cout << gas_species[i]->getSpeciesName() << std::endl;

  for (unsigned int i=0; i<gas_species.size(); ++i)
    gas_species[i]->init(spectral_indices);
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



