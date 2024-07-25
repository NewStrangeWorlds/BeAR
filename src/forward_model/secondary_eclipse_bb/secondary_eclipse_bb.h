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

#include "../../config/global_config.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../retrieval/priors.h"
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




class SecondaryEclipseBlackBodyModel : public ForwardModel{
  public:
    SecondaryEclipseBlackBodyModel (
      const SecondaryEclipseBlackBodyConfig model_config,
      Priors* priors_,
      GlobalConfig* config_,
      SpectralGrid* spectral_grid_,
      std::vector<Observation>& observations_);
    virtual ~SecondaryEclipseBlackBodyModel();
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
      const std::vector<double>& parameter,
      double* model_spectrum_gpu);
  protected:
    GlobalConfig* config;
    SpectralGrid* spectral_grid;

    StellarSpectrumModel* stellar_model;

    std::vector<Observation>& observations;
    size_t nb_observation_points = 0;

    size_t nb_general_param = 0;
    size_t nb_stellar_param = 0;
    size_t nb_spectrum_modifier_param = 0;

    size_t nb_total_param() {
      return nb_general_param 
             + nb_stellar_param
             + nb_spectrum_modifier_param;
    }

    virtual void setPriors(Priors* priors);
    void initModules(const SecondaryEclipseBlackBodyConfig& model_config);

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

    void postProcessSpectrum(
      std::vector<double>& model_spectrum, 
      std::vector<double>& model_spectrum_bands);
    void postProcessSpectrumGPU(
      double* model_spectrum, 
      double* model_spectrum_bands);

    void postProcessModel(
      const std::vector<double>& parameter,
      const std::vector<double>& model_spectrum_bands,
      std::vector<double>& temperature_profile,
      double& effective_temperature,
      std::vector<std::vector<double>>& mixing_ratios);

    bool testCPUvsGPU(const std::vector<double>& parameter, double* model_spectrum_gpu);
};


}


#endif

