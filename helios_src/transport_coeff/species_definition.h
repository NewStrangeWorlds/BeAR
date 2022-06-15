/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2022 Daniel Kitzmann
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

#include "opacity_species.h"

#include <vector>
#include <iostream>

#include "../chemistry/chem_species.h"
#include "../CUDA_kernels/cross_section_kernels.h"
#include "../additional/physical_const.h"
#include "../spectral_grid/spectral_grid.h"


namespace helios{


class GlobalConfig;
class SpectralGrid;


class GasGeneric : public OpacitySpecies {
  public:
    GasGeneric(
      GlobalConfig* config_ptr, 
      SpectralGrid* spectral_grid_ptr, 
      const unsigned int index, 
      const std::string name, 
      const std::string folder) 
        : OpacitySpecies(index, name, folder) 
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          init();
        }
    GasGeneric(
      GlobalConfig* config_ptr, 
      SpectralGrid* spectral_grid_ptr, 
      const unsigned int index, 
      const std::string name, 
      const std::string folder, 
      const size_t reference_species) 
        : OpacitySpecies(index, name, folder) 
        {
          config = config_ptr; spectral_grid = spectral_grid_ptr; 
          pressure_reference_species = reference_species; 
          init();
        }
    GasGeneric(
      GlobalConfig* config_ptr, 
      SpectralGrid* spectral_grid_ptr, 
      const unsigned int index, 
      const std::string name, 
      const std::string folder, 
      const std::vector<size_t> cia_collision_species) 
        : OpacitySpecies(index, name, folder)
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          cia_collision_partner = cia_collision_species; 
          init();
        }
    virtual ~GasGeneric() {}
};



class GasH : public OpacitySpecies {
  public:
    GasH(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr, const std::string folder) 
        : OpacitySpecies(_H, "H", folder)
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          init();
        }
    GasH(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr) 
        : OpacitySpecies(_H, "H", "")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          init();
        }
    virtual ~GasH() {}
  protected:
    virtual bool calcRalyleighCrossSections(std::vector<double>& cross_sections);
};



class GasHm : public OpacitySpecies {
  public:
    GasHm(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr) 
        : OpacitySpecies(_Hm, "H-", "")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          init();
        }
    virtual ~GasHm() {}
  protected:
    virtual void calcContinuumAbsorptionGPU(
      const double temperature, 
      const std::vector<double>& number_densities,
      const size_t nb_grid_points, 
      const size_t grid_point,
      double* absorption_coeff_device) 
      {
        calcHmContinuumHost(
          number_densities[_Hm],
          number_densities[_H],
          number_densities[_e]*constants::boltzmann_k*temperature,
          temperature,
          spectral_grid->nbSpectralPoints(), 
          grid_point,
          spectral_grid->wavelength_list_gpu,
          absorption_coeff_device);
      };
};



class GasH2 : public OpacitySpecies {
  public:
    GasH2(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr, const std::string folder) 
        : OpacitySpecies(_H2, "H2", folder)
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          init();
        }
    GasH2(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr) 
        : OpacitySpecies(_H2, "H2", "")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          init();
        }
    virtual ~GasH2() {}
  protected:
    virtual bool calcRalyleighCrossSections(std::vector<double>& cross_sections);
    virtual void calcRalyleighCrossSectionsGPU(
      const double number_density,
      const size_t nb_grid_points, 
      const size_t grid_point,
      double* scattering_coeff_dev);
};


class GasHe : public OpacitySpecies {
  public:
    GasHe(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr, const std::string folder)
        : OpacitySpecies(_He, "He", folder)
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          init();
        }
    GasHe(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr) 
        : OpacitySpecies(_He, "He", "")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          init();
        }
    virtual ~GasHe() {}
  protected:
    virtual bool calcRalyleighCrossSections(std::vector<double>& cross_sections);
    virtual void calcRalyleighCrossSectionsGPU(
      const double number_density,
      const size_t nb_grid_points, 
      const size_t grid_point,
      double* scattering_coeff_dev);
};

}

#endif
