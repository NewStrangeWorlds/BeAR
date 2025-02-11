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
#include "../generic_config.h"

#include "../../config/global_config.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../observations/observations.h"


namespace bear {


class FlatLinePostProcessConfig : public GenericConfig{
  public:
    bool save_spectra = false;
    bool delete_sampler_files = false;

    FlatLinePostProcessConfig (const std::string& folder_path);
    void readConfigFile(const std::string& file_name);
};



class FlatLine : public ForwardModel{
  public:
    FlatLine (
      GlobalConfig* config_, 
      SpectralGrid* spectral_grid_,
      std::vector<Observation>& observations_);
    virtual ~FlatLine();

    virtual size_t parametersNumber() {
      return nb_total_param();};
    
    virtual bool calcModelCPU(
      const std::vector<double>& parameter, 
      std::vector<double>& spectrum, 
      std::vector<std::vector<double>>& spectrum_obs);
    
    virtual bool calcModelGPU(
      const std::vector<double>& parameters, 
      double* spectrum, 
      std::vector<double*>& spectrum_obs);
    
    virtual void postProcess(
      const std::vector< std::vector<double> >& model_parameter,
      const size_t best_fit_model,
      bool& delete_unused_files);
    
    virtual bool testModel(
      const std::vector<double>& parameters);
  protected:
    size_t nb_general_param = 1;

    size_t nb_total_param() 
      {return nb_general_param;}

    void postProcessSpectrum(
      std::vector<double>& model_spectrum, std::vector<double>& model_spectrum_bands);
    void postProcessSpectrumGPU(
      double* model_spectrum, double* model_spectrum_bands);
};


}


#endif

