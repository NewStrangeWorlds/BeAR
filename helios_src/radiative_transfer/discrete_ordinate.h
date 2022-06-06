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


#ifndef _discrete_ordinate_h
#define _discrete_ordinate_h


#include <vector>
#include <iostream>
#include <cmath>

#include "radiative_transfer.h"
#include "../forward_model/atmosphere/atmosphere.h"


extern "C" {
  #include "../../cdisort_src/cdisort.h"
}


namespace helios {

//forward declaration
class Retrieval;



class DiscreteOrdinates : public RadiativeTransfer{
  public:
    DiscreteOrdinates(SpectralGrid* spectral_grid_ptr, const size_t nb_streams, const size_t nb_grid_points, const bool use_gpu);
    virtual ~DiscreteOrdinates() {finaliseDISORT();}
    virtual void calcSpectrum(const Atmosphere& atmosphere,
                              const std::vector< std::vector<double> >& absorption_coeff, 
                              const std::vector< std::vector<double> >& scattering_coeff,
                              const std::vector< std::vector<double> >& cloud_optical_depth,
                              const std::vector< std::vector<double> >& cloud_single_scattering,
                              const std::vector< std::vector<double> >& cloud_asym_param,
                              const double spectrum_scaling,
                              std::vector<double>& spectrum);
    virtual void calcSpectrumGPU(const Atmosphere& atmosphere,
                                 double* absorption_coeff_dev,
                                 double* scattering_coeff_dev,
                                 double* cloud_optical_depth,
                                 double* cloud_single_scattering,
                                 double* cloud_asym_param,
                                 const double spectrum_scaling,
                                 double* model_spectrum_dev) {std::cout << "Sorry, CDISORT has no GPU option :(\n";}
  private:
    SpectralGrid* spectral_grid;
    
    disort_state ds;
    disort_output out;

    double calcSpectrum(const std::vector<double> absorption_coeff, const std::vector<double> scattering_coeff, const std::vector<double>& cloud_optical_depth,
                        const std::vector<double>& vertical_grid, const size_t nu_index);
    void calcRadiativeTransfer(double incident_stellar_radiation, double zenith_angle,
                               std::vector<double>& flux_up, std::vector<double>& flux_down, std::vector<double>& mean_intensity);
    void receiveTemperatureStructure(const std::vector<double>& temperature_structure, const double& surface_temperature);
    void receiveTransportCoefficients(const double wavenumber_input, const std::vector<double>& optical_depth,
	                                                      const std::vector<double>& single_scattering_albedo, const std::vector<double>& asymmetry_parameter,
	                                                      const double surface_albedo);
    void initDISORT(unsigned int nb_streams, unsigned int nb_layers);
    void finaliseDISORT();
    void runDISORT(std::vector<double>& flux_up, std::vector<double>& flux_down, std::vector<double>& mean_intensity);
};


}
#endif


