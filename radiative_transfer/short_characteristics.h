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


#ifndef _short_char_h
#define _short_char_h


#include <vector>
#include <iostream>
#include <cmath>

#include "radiative_transfer.h"


namespace helios {

//forward declaration
class Retrieval;



class ShortCharacteristics : public RadiativeTransfer{
  public:
    ShortCharacteristics(Retrieval* retrieval_ptr) {retrieval = retrieval_ptr;}
    virtual ~ShortCharacteristics() {}
    virtual void calcSpectrum(const std::vector< std::vector<double> >& absorption_coeff, 
                              const std::vector< std::vector<double> >& scattering_coeff,
                              const std::vector<double>& temperature, 
                              const std::vector<double>& vertical_grid,
                              std::vector<double>& spectrum);
    virtual void calcSpectrumGPU(double* model_spectrum_dev,
                                 double* absorption_coeff_dev, 
                                 double* scattering_coeff_dev, 
                                 double* wavenumber_list_dev,
                                 const std::vector<double>& cloud_optical_depth,
                                 const std::vector<double>& temperature, 
                                 const std::vector<double>& vertical_grid,
                                 const double radius_distance_scaling);

  private:
    Retrieval* retrieval;

    const std::vector<double> gauss_nodes{0.339981, 0.861136};
    const std::vector<double> gauss_weights{0.652145, 0.347855};
    const size_t nb_angles = gauss_nodes.size();

    double calcSpectrum(const std::vector<double>& absorption_coeff, const std::vector<double>& temperature, const std::vector<double>& vertical_grid,
                        const size_t nu_index);

};


}




#endif


