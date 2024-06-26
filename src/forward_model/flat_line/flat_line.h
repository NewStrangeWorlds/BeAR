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


namespace helios {


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
      const std::vector< std::vector<double> >& model_spectrum_bands,
      const size_t best_fit_model);

    virtual std::vector<double> convertSpectrumToModel(const std::vector<double>& spectrum);
    
    virtual bool testModel(
      const std::vector<double>& parameter, double* model_spectrum_gpu);
  protected:
    GlobalConfig* config;
    SpectralGrid* spectral_grid;

    std::vector<Observation>& observations;
    size_t nb_observation_points = 0;
    
    size_t nb_general_param = 1;

    size_t nb_total_param() 
      {return nb_general_param;}

    virtual void setPriors(Priors* priors);
    void readPriorConfigFile(
      const std::string& file_name, 
      std::vector<std::string>& prior_type, 
      std::vector<std::string>& prior_description, 
      std::vector<std::vector<double>>& prior_parameter);

    void postProcessSpectrum(
      std::vector<double>& model_spectrum, std::vector<double>& model_spectrum_bands);
    void postProcessSpectrumGPU(
      double* model_spectrum, double* model_spectrum_bands);
    bool testCPUvsGPU(
      const std::vector<double>& parameter, double* model_spectrum_gpu);
};


}


#endif

