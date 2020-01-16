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


#ifndef _brown_dward_h
#define _brown_dwarf_h

#include <vector>
#include <iostream>
#include <cmath>

#include "forward_model.h"
#include "../transport_coeff/transport_coeff.h"
#include "../radiative_transfer/discrete_ordinate.h"
#include "../radiative_transfer/short_characteristics.h"
#include "piecewise_poly.h"

#include "../radiative_transfer/radiative_transfer.h"


namespace helios {


class Retrieval;



class BrownDwarfModel : public ForwardModel{
  public:
    BrownDwarfModel (Retrieval* retrieval_ptr, const size_t& nb_points, const double domain_boundaries [2]);
    virtual ~BrownDwarfModel();
    bool calcModel(const std::vector<double>& parameter, std::vector<double>& spectrum);
    bool calcModelGPU(const std::vector<double>& parameter, double* model_spectrum);

    std::vector<double> temperatureProfile(const std::vector<double>& parameter, std::vector<double>& pressure_profile);
    double radiusDistanceScaling(const std::vector<double>& parameter);

  protected:
    TransportCoefficients transport_coeff;
    RadiativeTransfer* radiative_transfer = nullptr;
    
    size_t nb_grid_points = 0;

    std::vector<double> pressure;
    std::vector<double> temperature;
    std::vector<double> z_grid;
    std::vector< std::vector<double> > number_densities;

    double radius_distance_scaling = 0;

    PiecewisePolynomial temperature_pol;

    std::vector< std::vector<double> > absorption_coeff;
    std::vector< std::vector<double> > scattering_coeff;
    std::vector<double> cloud_optical_depths;

    //pointer to the array that holds the pointers to the coefficients on the GPU
    double* absorption_coeff_gpu = nullptr;

    void createPressureGrid(const double domain_boundaries [2]);
    bool calcAtmosphereStructure(const std::vector<double>& parameter);

    void calcTemperature (const std::vector<double>& temp_parameters);
    void calcFreeChemistry (const std::vector<double>& chem_parameters, std::vector<double>& mean_molecular_weight);
    void calcAltitude(const double g, const std::vector<double>& mean_molecular_weights);

    void calcGreyCloudLayer(const std::vector<double>& cloud_parameters);
    void calcCloudPosition(const double top_pressure, const double bottom_pressure, unsigned int& top_index, unsigned int& bottom_index);
};


}


#endif

