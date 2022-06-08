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


#ifndef TRANSPORT_COEFF_H
#define TRANSPORT_COEFF_H

#include "opacity_species.h"


#include <vector>



namespace helios{


//forward declaration
class GlobalConfig;



class TransportCoefficients {
  public:
    TransportCoefficients(GlobalConfig* config_ptr, SpectralGrid* grid_ptr, 
                          const std::vector<std::string>& opacity_species_symbol, const std::vector<std::string>& opacity_species_folder);
    ~TransportCoefficients();

    void calcTransportCoefficients(const double temperature, const double pressure, const std::vector<double>& number_densities,
                                   std::vector<double>& absorption_coeff, std::vector<double>& scattering_coeff);

    void calcTransportCoefficientsGPU(const double temperature, const double pressure, const std::vector<double>& number_densities,
                                      const size_t nb_grid_points, const size_t grid_point,
                                      double* absorption_coeff_device, double* scattering_coeff_device);
    
  private:
    GlobalConfig* config = nullptr;
    SpectralGrid* spectral_grid = nullptr;

    std::vector<OpacitySpecies*> gas_species;

    bool addOpacitySpecies(const std::string& species_symbol, const std::string& species_folder);
};


}

#endif
