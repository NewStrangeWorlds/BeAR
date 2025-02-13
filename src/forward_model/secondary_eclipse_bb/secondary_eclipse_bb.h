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


#ifndef _secondary_eclipse_bb_h
#define _secondary_eclipse_bb_h

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

#include "../stellar_spectrum/stellar_spectrum.h"


namespace bear {


//this struct handles config
//it will read in the corresponding parameter file
//and will then be used to create a model object
struct SecondaryEclipseBlackBodyConfig{
  std::string stellar_spectrum_model;
  std::vector<std::string> stellar_model_parameters;

  SecondaryEclipseBlackBodyConfig (const std::string& folder_path);
  void readConfigFile(const std::string& file_name);
};



class SecondaryEclipseBlackBodyPostConfig : public GenericConfig{
  public:
    bool save_spectra = true;
    bool delete_sampler_files = false;

    SecondaryEclipseBlackBodyPostConfig (const std::string& folder_path);
    void readConfigFile(const std::string& file_name);
};



class SecondaryEclipseBlackBodyModel : public ForwardModel{
  public:
    SecondaryEclipseBlackBodyModel (
      const SecondaryEclipseBlackBodyConfig model_config,
      GlobalConfig* config_,
      SpectralGrid* spectral_grid_,
      std::vector<Observation>& observations_);
    virtual ~SecondaryEclipseBlackBodyModel();
    
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
    virtual void postProcess(
      GenericConfig* post_process_config_,
      const std::vector< std::vector<double> >& model_parameter,
      const size_t best_fit_model,
      bool& delete_unused_files) {};

    virtual bool testModel(
      const std::vector<double>& parameters);
  protected:
    StellarSpectrumModel* stellar_model;

    size_t nb_general_param = 0;
    size_t nb_stellar_param = 0;

    size_t nb_total_param() {
      return nb_general_param 
             + nb_stellar_param
             + nb_spectrum_modifier_param;
    }

    void initModules(const SecondaryEclipseBlackBodyConfig& model_config);

    std::vector<double> model_parameters;
    std::vector<double> stellar_parameters;
    std::vector<double> spectrum_modifier_parameters;

    void extractParameters(
      const std::vector<double>& parameters);

    void calcSecondaryEclipseGPU(
      double* secondary_eclipse,
      double* planet_spectrum,
      const double* stellar_spectrum,
      const int nb_points,
      const double radius_ratio,
      const double* albedo_contribution);
    
    void calcPlanetSpectrumGPU(
      const double planet_temperature,
      double* spectrum_dev);

    void postProcessModel(
      const std::vector<double>& parameter,
      const std::vector<double>& model_spectrum_bands,
      std::vector<double>& temperature_profile,
      double& effective_temperature,
      std::vector<std::vector<double>>& mixing_ratios);
};


}


#endif

