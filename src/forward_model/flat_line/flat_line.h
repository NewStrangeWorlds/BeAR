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


#ifndef _flat_line_h
#define _flat_line_h

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

#include "../forward_model.h"

#include "../../config/global_config.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../retrieval/priors.h"
#include "../../observations/observations.h"


namespace bear {


class FlatLine : public ForwardModel{
  public:
    FlatLine (
      Priors* priors_,
      GlobalConfig* config_, 
      SpectralGrid* spectral_grid_,
      std::vector<Observation>& observations_);
    virtual ~FlatLine();
    
    virtual bool calcModel(
      const std::vector<double>& parameter, 
      std::vector<double>& spectrum, 
      std::vector<double>& model_spectrum_bands);
    
    virtual bool calcModelGPU(
      const std::vector<double>& parameter, 
      double* model_spectrum, 
      double* model_spectrum_bands);
    
    virtual void postProcess(
      const std::vector< std::vector<double> >& model_parameter,
      const size_t best_fit_model);
    
    virtual bool testModel(
      const std::vector<double>& parameter, double* model_spectrum_gpu);
  protected:
    size_t nb_general_param = 1;

    size_t nb_total_param() 
      {return nb_general_param;}

    virtual void setPriors(Priors* priors);

    void postProcessSpectrum(
      std::vector<double>& model_spectrum, std::vector<double>& model_spectrum_bands);
    void postProcessSpectrumGPU(
      double* model_spectrum, double* model_spectrum_bands);
    bool testCPUvsGPU(
      const std::vector<double>& parameter, double* model_spectrum_gpu);
};


}


#endif

