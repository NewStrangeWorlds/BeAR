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


#ifndef SPECIES_DEFINITION_H
#define SPECIES_DEFINITION_H

#include "transport_coeff_single_species.h"

#include <vector>
#include <iostream>

#include "../chemistry/chem_species.h"

namespace helios{


class GlobalConfig;
class SpectralGrid;


//generic class for an opacity species that does not require any special treatment
//i.e. no additional continuum absorption or Rayleigh scattering coefficients
//this class will just handle the line absorption
class GasGeneric : public TransportCoefficientsSingleSpecies {
  public:
    GasGeneric(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr, const unsigned int index, const std::string name) 
      {config = config_ptr; spectral_grid = spectral_grid_ptr; species_index = index; species_name = name;}
    virtual ~GasGeneric() {}
  protected:
    virtual bool calcContinuumAbsorption(const double temperature, const std::vector<double>& number_densities, std::vector<double>& absorption_coeff) {return false;}
    virtual void calcContinuumAbsorptionGPU(const double temperature, const std::vector<double>& number_densities,
                                            const size_t nb_grid_points, const size_t grid_point,
                                            double* absorption_coeff_device) {}
    virtual void calcRalyleighCrossSections(std::vector<double>& cross_sections) {}
};



//separate class for H2
//contains the HITRAN H2-H2 and H2-H CIA
class GasH2 : public TransportCoefficientsSingleSpecies {
  public:
    GasH2(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr) {config = config_ptr; spectral_grid = spectral_grid_ptr; species_index = _H2; species_name = "H2";}
    virtual ~GasH2() {}
  protected:
    virtual bool calcContinuumAbsorption(const double temperature, const std::vector<double>& number_densities, std::vector<double>& absorption_coeff);
    virtual void calcContinuumAbsorptionGPU(const double temperature, const std::vector<double>& number_densities,
                                            const size_t nb_grid_points, const size_t grid_point,
                                            double* absorption_coeff_device);
    virtual void calcRalyleighCrossSections(std::vector<double>& cross_sections);
  private:
    CIACoefficients cia_H2H2;
    CIACoefficients cia_H2He;

    void readHITRAN_H2H2_CIA();
    void readHITRAN_H2He_CIA();
};



}

#endif
