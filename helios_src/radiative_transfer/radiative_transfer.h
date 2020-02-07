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


#ifndef _radiative_transfer_h
#define _radiative_transfer_h

#include <vector>



namespace helios {

//forward declaration
class SpectralGrid;



class RadiativeTransfer{
  public:
    virtual ~RadiativeTransfer() {}
    virtual void calcSpectrum(const std::vector< std::vector<double> >& absorption_coeff, 
                              const std::vector< std::vector<double> >& scattering_coeff,
                              const std::vector<double>& cloud_optical_depth,
                              const std::vector<double>& temperature, 
                              const std::vector<double>& vertical_grid,
                              std::vector<double>& spectrum) = 0;
    virtual void calcSpectrumGPU(double* model_spectrum_dev,
                                 double* absorption_coeff_dev, 
                                 double* scattering_coeff_dev, 
                                 double* wavenumber_list_dev,
                                 const std::vector<double>& cloud_optical_depth,
                                 const std::vector<double>& temperature, 
                                 const std::vector<double>& vertical_grid,
                                 const double radius_distance_scaling) = 0;
};


}
#endif

