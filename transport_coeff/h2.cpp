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


#include "species_definition.h"
#include "transport_coeff_single_species.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <assert.h>
#include <omp.h>


#include "../config/global_config.h"
#include "../spectral_grid/spectral_grid.h"





namespace helios{


void GasH2::readHITRAN_H2H2_CIA()
{
  std::string file_path = config->crossSectionFilePath();
  std::string file_name = file_path + species_name+"/H2-H2_2011.cia";


  std::cout << "reading H2-H2 CIA file: " << file_name << std::endl;
  cia_H2H2.init(file_name, spectral_grid->wavenumber_list, config->useGPU());


  cia_H2H2.species_id_1 = _H2;
  cia_H2H2.species_id_2 = _H2;

}



void GasH2::readHITRAN_H2He_CIA()
{
  std::string file_path = config->crossSectionFilePath();
  std::string file_name = file_path + species_name+"/H2-He_2011a.cia";


  std::cout << "reading H2-He CIA file: " << file_name << std::endl;
  cia_H2He.init(file_name, spectral_grid->wavenumber_list, config->useGPU());


  cia_H2He.species_id_1 = _H2;
  cia_H2He.species_id_2 = _He;

}




bool GasH2::calcContinuumAbsorption(const double temperature, const std::vector<double>& number_densities, std::vector<double>& absorption_coeff)
{

  if (cia_H2H2.nb_temperatures == 0)
    readHITRAN_H2H2_CIA();

  if (cia_H2He.nb_temperatures == 0)
    readHITRAN_H2He_CIA();


  std::vector<double> log_cross_section;


  if (number_densities[_H2] > 0)
  {
    log_cross_section = cia_H2H2.calcCIACoefficients(temperature);

    for (unsigned int j=0; j<log_cross_section.size(); j++)
      absorption_coeff[j] += std::pow(10.0, log_cross_section[j]) * number_densities[_H2] * number_densities[_H2];
  }



  if (number_densities[_H2] > 0 && number_densities[_He] > 0)
  {
    log_cross_section = cia_H2He.calcCIACoefficients(temperature);

    for (unsigned int j=0; j<log_cross_section.size(); j++)
      absorption_coeff[j] += std::pow(10., log_cross_section[j]) * number_densities[_H2] * number_densities[_He];
  }



  return true;
}




void GasH2::calcContinuumAbsorptionGPU(const double temperature, const std::vector<double>& number_densities,
                                       const size_t nb_grid_points, const size_t grid_point,
                                       double* absorption_coeff_device)
{

  if (cia_H2H2.nb_temperatures == 0)
    readHITRAN_H2H2_CIA();

  if (cia_H2He.nb_temperatures == 0)
    readHITRAN_H2He_CIA();


  if (number_densities[_H2] > 0)
    cia_H2H2.calcCIACoefficientsGPU(temperature, number_densities[_H2]*number_densities[_H2],
                                    nb_grid_points, grid_point,
                                    absorption_coeff_device);


  if (number_densities[_H2] > 0 && number_densities[_He] > 0)
    cia_H2He.calcCIACoefficientsGPU(temperature, number_densities[_H2]*number_densities[_He],
                                    nb_grid_points, grid_point,
                                    absorption_coeff_device);

}




void GasH2::calcRalyleighCrossSections(std::vector<double>& cross_sections)
{



}



}

