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


#ifndef SPECIES_DEFINITION_H
#define SPECIES_DEFINITION_H

#include "opacity_species.h"

#include <vector>
#include <iostream>

#include "../chemistry/chem_species.h"
#include "../CUDA_kernels/cross_section_kernels.h"
#include "../additional/physical_const.h"
#include "../spectral_grid/spectral_grid.h"


namespace bear{


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


class GasHm : public OpacitySpecies {
  public:
    GasHm(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr) 
        : OpacitySpecies(_Hm, "H-", "Continuum")
        { 
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          continuum_available = true;
          init();
        }
    virtual ~GasHm() {}
  protected:
    virtual bool calcContinuumAbsorption(
      const double temperature,
      const std::vector<double>& number_densities,
      std::vector<double>& absorption_coeff);
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
    private:
      std::vector<double> boundFreeAbsorption(const double temperature);
      std::vector<double> freeFreeAbsorption(const double temperature);
};


class GasHRayleigh : public OpacitySpecies {
  public:
    GasHRayleigh(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr, const std::string folder) 
        : OpacitySpecies(_H, "H Rayleigh", "Rayleigh")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          rayleigh_available = true;
          init();
        }
    GasHRayleigh(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr) 
        : OpacitySpecies(_H, "H Rayleigh", "Rayleigh")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          rayleigh_available = true;
          init();
        }
    virtual ~GasHRayleigh() {}
  protected:
    virtual bool calcRayleighCrossSections(std::vector<double>& cross_sections);
    virtual void calcRayleighCrossSectionsGPU(
      const double number_density,
      const size_t nb_grid_points, 
      const size_t grid_point,
      double* scattering_coeff_dev);
};


class GasH2Rayleigh : public OpacitySpecies {
  public:
    GasH2Rayleigh(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr, const std::string folder) 
        : OpacitySpecies(_H2, "H2 Rayleigh", "Rayleigh")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          rayleigh_available = true;
          init();
        }
    GasH2Rayleigh(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr) 
        : OpacitySpecies(_H2, "H2 Rayleigh", "Rayleigh")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          rayleigh_available = true;
          init();
        }
    virtual ~GasH2Rayleigh() {}
  protected:
    virtual bool calcRayleighCrossSections(std::vector<double>& cross_sections);
    virtual void calcRayleighCrossSectionsGPU(
      const double number_density,
      const size_t nb_grid_points, 
      const size_t grid_point,
      double* scattering_coeff_dev);
};


class GasHeRayleigh : public OpacitySpecies {
  public:
    GasHeRayleigh(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr, const std::string folder)
        : OpacitySpecies(_He, "He Rayleigh", "Rayleigh")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          rayleigh_available = true;
          init();
        }
    GasHeRayleigh(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr) 
        : OpacitySpecies(_He, "He Rayleigh", "Rayleigh")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          rayleigh_available = true;
          init();
        }
    virtual ~GasHeRayleigh() {}
  protected:
    virtual bool calcRayleighCrossSections(std::vector<double>& cross_sections);
    virtual void calcRayleighCrossSectionsGPU(
      const double number_density,
      const size_t nb_grid_points, 
      const size_t grid_point,
      double* scattering_coeff_dev);
};


class GasCORayleigh : public OpacitySpecies {
  public:
    GasCORayleigh(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr, const std::string folder)
        : OpacitySpecies(_CO, "CO Rayleigh", "Rayleigh")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          rayleigh_available = true;
          init();
        }
    GasCORayleigh(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr) 
        : OpacitySpecies(_CO, "CO Rayleigh", "Rayleigh")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          rayleigh_available = true;
          init();
        }
    virtual ~GasCORayleigh() {}
  protected:
    virtual bool calcRayleighCrossSections(std::vector<double>& cross_sections);
    virtual void calcRayleighCrossSectionsGPU(
      const double number_density,
      const size_t nb_grid_points, 
      const size_t grid_point,
      double* scattering_coeff_dev);
};


class GasCO2Rayleigh : public OpacitySpecies {
  public:
    GasCO2Rayleigh(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr, const std::string folder)
        : OpacitySpecies(_CO2, "CO2 Rayleigh", "Rayleigh")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          rayleigh_available = true;
          init();
        }
    GasCO2Rayleigh(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr) 
        : OpacitySpecies(_CO2, "CO2 Rayleigh", "Rayleigh")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          rayleigh_available = true;
          init();
        }
    virtual ~GasCO2Rayleigh() {}
  protected:
    virtual bool calcRayleighCrossSections(std::vector<double>& cross_sections);
    virtual void calcRayleighCrossSectionsGPU(
      const double number_density,
      const size_t nb_grid_points, 
      const size_t grid_point,
      double* scattering_coeff_dev);
};


class GasCH4Rayleigh : public OpacitySpecies {
  public:
    GasCH4Rayleigh(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr, const std::string folder)
        : OpacitySpecies(_CH4, "CH4 Rayleigh", "Rayleigh")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          rayleigh_available = true;
          init();
        }
    GasCH4Rayleigh(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr) 
        : OpacitySpecies(_CH4, "CH4 Rayleigh", "Rayleigh")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          rayleigh_available = true;
          init();
        }
    virtual ~GasCH4Rayleigh() {}
  protected:
    virtual bool calcRayleighCrossSections(std::vector<double>& cross_sections);
    virtual void calcRayleighCrossSectionsGPU(
      const double number_density,
      const size_t nb_grid_points, 
      const size_t grid_point,
      double* scattering_coeff_dev);
};



class GasH2ORayleigh : public OpacitySpecies {
  public:
    GasH2ORayleigh(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr, const std::string folder)
        : OpacitySpecies(_H2O, "H2O Rayleigh", "Rayleigh")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          rayleigh_available = true;
          init();
        }
    GasH2ORayleigh(GlobalConfig* config_ptr, SpectralGrid* spectral_grid_ptr) 
        : OpacitySpecies(_H2O, "H2O Rayleigh", "Rayleigh")
        {
          config = config_ptr; 
          spectral_grid = spectral_grid_ptr; 
          rayleigh_available = true;
          init();
        }
    virtual ~GasH2ORayleigh() {}
  protected:
    virtual bool calcRayleighCrossSections(std::vector<double>& cross_sections);
    virtual void calcRayleighCrossSectionsGPU(
      const double number_density,
      const size_t nb_grid_points, 
      const size_t grid_point,
      double* scattering_coeff_dev);
};





}

#endif
