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


#ifndef _emission_h
#define _emission_h

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
#include "../atmosphere/atmosphere.h"
#include "../../cloud_model/cloud_model.h"
#include "../../chemistry/chemistry.h"
#include "../../temperature/temperature.h"
#include "../../transport_coeff/transport_coeff.h"
#include "../../radiative_transfer/radiative_transfer.h"
#include "../../transport_coeff/opacity_calc.h"


namespace bear {


//this struct handles the Emission config
//it will read in the corresponding parameter file
//and will then be used to create a model object
struct EmissionModelConfig{
  size_t nb_grid_points = 0;

  double atmos_boundaries[2] {0, 0};
  double atmos_top_pressure = 0;
  double atmos_bottom_pressure = 0;

  bool use_cloud_model = false;

  std::string temperature_profile_model;
  std::vector<std::string> temperature_profile_parameters;

  std::string radiative_transfer_model;
  std::vector<std::string> radiative_transfer_parameters;

  std::vector<std::string> chemistry_model;
  std::vector<std::vector<std::string>> chemistry_parameters;

  std::vector<std::string> cloud_model;
  std::vector<std::vector<std::string>> cloud_model_parameters;

  std::vector<std::string> opacity_species_symbol;
  std::vector<std::string> opacity_species_folder;

  EmissionModelConfig (const std::string& folder_path);
  void readConfigFile(const std::string& file_name);
  void readCloudConfig(std::fstream& file);
  void readChemistryConfig(std::fstream& file);
  void readOpacityConfig(std::fstream& file);
};


class EmissionModel : public ForwardModel{
  public:
    EmissionModel (
      const EmissionModelConfig model_config,
      Priors* priors_,
      GlobalConfig* config_, 
      SpectralGrid* spectral_grid_,
      std::vector<Observation>& observations_);
    virtual ~EmissionModel();
    
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

    Atmosphere atmosphere;
    OpacityCalculation opacity_calc;

    RadiativeTransfer* radiative_transfer = nullptr;
    Temperature* temperature_profile = nullptr;
    std::vector<Chemistry*> chemistry;
    std::vector<CloudModel*> cloud_models;

    std::vector<Observation>& observations;
    size_t nb_observation_points = 0;
    
    size_t nb_general_param = 0;
    size_t nb_total_chemistry_param = 0;
    size_t nb_temperature_param = 0;
    size_t nb_total_cloud_param = 0;
    size_t nb_spectrum_modifier_param = 0;

    size_t nb_total_param() 
      {return nb_general_param 
        + nb_total_chemistry_param 
        + nb_temperature_param 
        + nb_total_cloud_param
        + nb_spectrum_modifier_param;}
    
    size_t nb_grid_points = 0;

    virtual void setPriors(Priors* priors);
    void initModules(const EmissionModelConfig& model_config);

    bool calcAtmosphereStructure(const std::vector<double>& parameter);
    double radiusDistanceScaling(const std::vector<double>& parameter);

    void postProcessSpectrum(
      std::vector<double>& model_spectrum, std::vector<double>& model_spectrum_bands);
    void postProcessSpectrumGPU(
      double* model_spectrum, double* model_spectrum_bands);

    void postProcessModel(
      const std::vector<double>& parameter, 
      const std::vector<double>& model_spectrum_bands, 
      std::vector<double>& temperature_profile, 
      double& effective_temperature,
      std::vector<std::vector<double>>& mixing_ratios);
    double postProcessEffectiveTemperature(
      const std::vector<double>& model_spectrum_bands, const double radius_distance_scaling);
    void savePostProcessChemistry(
      const std::vector<std::vector<std::vector<double>>>& mixing_ratios, 
      const unsigned int species);
    void savePostProcessTemperatures(
      const std::vector<std::vector<double>>& temperature_profiles);
    void savePostProcessEffectiveTemperatures(
      const std::vector<double>& effective_temperatures);
    void postProcessContributionFunctions(
      const std::vector<double>& model_parameter);
    void saveContributionFunctions(
      std::vector< std::vector<double>>& contribution_function,
      const size_t observation_index);

    bool testCPUvsGPU(
      const std::vector<double>& parameter, double* model_spectrum_gpu);
};


}


#endif

