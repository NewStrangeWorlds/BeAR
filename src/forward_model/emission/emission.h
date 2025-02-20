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
#include "../generic_config.h"

#include "../../config/global_config.h"
#include "../../spectral_grid/spectral_grid.h"
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
class EmissionModelConfig : public GenericConfig {
  public:
    size_t nb_grid_points = 0;

    std::vector<double> atmos_boundaries = {0, 0};

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

    EmissionModelConfig (
      const std::string& folder_path);
    EmissionModelConfig (
      const int nb_grid_points_,
      const double atmos_bottom_pressure_,
      const double atmos_top_pressure_,
      const std::string& temperature_profile_model_,
      const std::vector<std::string>& temperature_profile_parameters_,
      const std::string radiative_transfer_model_,
      const std::vector<std::string>& radiative_transfer_parameters_,
      const std::vector<std::string>& chemistry_model_,
      const std::vector<std::vector<std::string>>& chemistry_parameters_,
      const std::vector<std::string>& opacity_species_symbol_,
      const std::vector<std::string>& opacity_species_folder_);
    EmissionModelConfig (
      const int nb_grid_points_,
      const double atmos_bottom_pressure_,
      const double atmos_top_pressure_,
      const std::string& temperature_profile_model_,
      const std::vector<std::string>& temperature_profile_parameters_,
      const std::string radiative_transfer_model_,
      const std::vector<std::string>& radiative_transfer_parameters_,
      const std::vector<std::string>& chemistry_model_,
      const std::vector<std::vector<std::string>>& chemistry_parameters_,
      const std::vector<std::string>& opacity_species_symbol_,
      const std::vector<std::string>& opacity_species_folder_,
      const std::vector<std::string>& cloud_model_,
      const std::vector<std::vector<std::string>>& cloud_model_parameters_);
    
    void readConfigFile(const std::string& file_name);
};


class EmissionPostProcessConfig : public GenericConfig{
  public:
    bool save_temperatures = true;
    bool save_effective_temperatures = true;
    bool save_spectra = true;
    bool save_contribution_functions = false;
    bool delete_sampler_files = false;
    std::vector<chemical_species_id> species_to_save;

    EmissionPostProcessConfig (
      const std::string& folder_path);
    EmissionPostProcessConfig (
      const bool save_temperatures_, 
      const bool save_effective_temperatures_, 
      const bool save_spectra_, 
      const bool save_contribution_functions_,
      const std::vector<std::string>& species_to_save_);
    
    void readConfigFile(const std::string& file_name);
};


class EmissionModel : public ForwardModel{
  public:
    EmissionModel (
      const EmissionModelConfig model_config,
      GlobalConfig* config_, 
      SpectralGrid* spectral_grid_,
      std::vector<Observation>& observations_);
    EmissionModel (
      GlobalConfig* config_, 
      SpectralGrid* spectral_grid_,
      const size_t nb_grid_points_,
      const std::string radiative_transfer_desc,
      const std::vector<std::string>& radiative_transfer_param,
      const std::vector<std::string>& opacity_species_symbol,
      const std::vector<std::string>& opacity_species_folder);
    virtual ~EmissionModel();

    virtual size_t parametersNumber() {
      return nb_total_param();};
    
    virtual bool calcModelCPU(
      const std::vector<double>& parameter, 
      std::vector<double>& spectrum, 
      std::vector<std::vector<double>>& spectrum_obs);
    
    virtual bool calcModelGPU(
      const std::vector<double>& parameter, 
      double* spectrum, 
      std::vector<double*>& spectrum_obs);
    
    std::vector<double> calcSpectrum(
      const double surface_gravity,
      const double radius,
      const double distance,
      const std::vector<double>& pressure,
      const std::vector<double>& temperature,
      const std::vector<std::string>& species_symbol,
      const std::vector<std::vector<double>>& mixing_ratios,
      const std::vector<std::vector<double>>& cloud_optical_depth,
      const bool use_variable_gravity);

    virtual void postProcess(
      const std::vector< std::vector<double> >& model_parameter,
      const size_t best_fit_model,
      bool& delete_unused_files);
    virtual void postProcess(
      GenericConfig* post_process_config_,
      const std::vector< std::vector<double> >& model_parameter,
      const size_t best_fit_model,
      bool& delete_unused_files);
    
    virtual bool testModel(
      const std::vector<double>& parameters);
  protected:
    Atmosphere atmosphere;
    OpacityCalculation opacity_calc;

    RadiativeTransfer* radiative_transfer = nullptr;
    Temperature* temperature_profile = nullptr;
    std::vector<Chemistry*> chemistry;
    std::vector<CloudModel*> cloud_models;
    
    size_t nb_general_param = 0;
    size_t nb_total_chemistry_param = 0;
    size_t nb_temperature_param = 0;
    size_t nb_total_cloud_param = 0;
   
    size_t nb_total_param() 
      {return nb_general_param 
        + nb_total_chemistry_param 
        + nb_temperature_param 
        + nb_total_cloud_param
        + nb_spectrum_modifier_param;}
    
    size_t nb_grid_points = 0;

    void initModules(const EmissionModelConfig& model_config);

    std::vector<double> model_parameters;
    std::vector<double> chemistry_parameters;
    std::vector<double> cloud_parameters;
    std::vector<double> temperature_parameters;
    std::vector<double> spectrum_modifier_parameters;

    void extractParameters(
      const std::vector<double>& parameters);

    bool calcAtmosphereStructure(const std::vector<double>& parameter);
    double radiusDistanceScaling(const std::vector<double>& parameter);
    void changeSpectrumUnitsGPU(double* spectrum_gpu);

    void setCloudProperties(const std::vector<std::vector<double>>& cloud_optical_depth);

    void postProcess(
      const EmissionPostProcessConfig& post_process_config_,
      const std::vector< std::vector<double> >& model_parameter,
      const size_t best_fit_model);
    void postProcessModel(
      const std::vector<double>& parameter, 
      const double integrated_flux, 
      std::vector<double>& temperature_profile, 
      double& effective_temperature,
      std::vector<std::vector<double>>& mixing_ratios);
    void calcPostProcessSpectra(
      const std::vector< std::vector<double> >& model_parameter,
      const size_t best_fit_model,
      const bool save_spectra,
      std::vector<double>& integrated_flux);
    double postProcessEffectiveTemperature(
      const double integrated_flux, 
      const double radius_distance_scaling);
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
};


}


#endif

