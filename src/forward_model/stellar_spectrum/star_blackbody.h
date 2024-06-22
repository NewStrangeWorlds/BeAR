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


#ifndef _star_blackbody_h
#define _star_blackbody_h

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

#include "stellar_spectrum.h"
#include "../../spectral_grid/spectral_grid.h"


namespace bear {


class StarBlackBody : public StellarSpectrumModel{
  public:
    StarBlackBody (SpectralGrid* spectral_grid_) 
      : spectral_grid(spectral_grid_)
      {nb_parameters = 1;}
    virtual ~StarBlackBody() {}
    
    virtual std::vector<double> calcFlux(
      const std::vector<double>& parameter);
    
    virtual void calcFluxGPU(
      const std::vector<double>& parameter,
      double* spectrum_gpu);
  protected:
    SpectralGrid* spectral_grid;
};


}

#endif