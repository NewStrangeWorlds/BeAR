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


#ifndef _transmission_h
#define _transmission_h


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
#include "../../transport_coeff/opacity_calc.h"
#include "../modules/module.h"
#include "../../chemistry/chem_species.h"


namespace bear {


constexpr double transmission_optical_depth_cutoff = 10;
constexpr double transmission_cutoff = 4.5399929e-5; //exp(-tau)


//this struct handles the Transmission Spectrum config
//it will read in the corresponding parameter file
//and will then be used to create a model object
class TransmissionModelConfig : public GenericConfig{
  public:
    size_t nb_grid_points = 0;

    std::vector<double> atmos_boundaries = {0, 0};

    bool fit_mean_molecular_weight = false;
    bool fit_scale_height = false;
    bool use_variable_gravity = false;

    //bool use_optional_modules = false;

    std::string temperature_profile_model;
    std::vector<std::string> temperature_profile_parameters;

    std::vector<std::string> chemistry_model;
    std::vector<std::vector<std::string>> chemistry_parameters;

    std::vector<std::string> cloud_model;
    std::vector<std::vector<std::string>> cloud_model_parameters;

    std::vector<std::string> modules;
    std::vector<std::vector<std::string>> modules_parameters;

    std::vector<std::string> opacity_species_symbol;
    std::vector<std::string> opacity_species_folder;
    
    TransmissionModelConfig (const std::string& folder_path);
    TransmissionModelConfig (
      const int nb_grid_points_,
      const double atmos_bottom_pressure_,
      const double atmos_top_pressure_,
      const std::string& temperature_profile_model_,
      const std::vector<std::string>& temperature_profile_parameters_,
      const std::vector<std::string>& chemistry_model_,
      const std::vector<std::vector<std::string>>& chemistry_parameters_,
      const std::vector<std::string>& opacity_species_symbol_,
      const std::vector<std::string>& opacity_species_folder_);
    TransmissionModelConfig (
      const bool fit_mean_molecular_weight_, 
      const bool fit_scale_height_, 
      const bool use_variable_gravity_,
      const int nb_grid_points_,
      const double atmos_bottom_pressure_,
      const double atmos_top_pressure_,
      const std::string& temperature_profile_model_,
      const std::vector<std::string>& temperature_profile_parameters_,
      const std::vector<std::string>& chemistry_model_,
      const std::vector<std::vector<std::string>>& chemistry_parameters_,
      const std::vector<std::string>& opacity_species_symbol_,
      const std::vector<std::string>& opacity_species_folder_,
      const std::vector<std::string>& cloud_model_,
      const std::vector<std::vector<std::string>>& cloud_model_parameters_,
      const std::vector<std::string>& modules_,
      const std::vector<std::vector<std::string>>& modules_parameters_);
    virtual void readConfigFile(const std::string& file_name);
};


class TransmissionPostProcessConfig : public GenericConfig{
  public:
    bool save_temperatures = false;
    bool save_spectra = true;
    bool delete_sampler_files = false;
    std::vector<chemical_species_id> species_to_save;
    
    TransmissionPostProcessConfig (
      const std::string& folder_path);
    TransmissionPostProcessConfig (
      const bool save_temperatures_, 
      const bool save_spectra_, 
      const std::vector<std::string>& species_to_save_);
    void readConfigFile(const std::string& file_name);
};



class TransmissionModel : public ForwardModel{
  public:
    TransmissionModel (
      const TransmissionModelConfig model_config,
      GlobalConfig* config_, 
      SpectralGrid* spectral_grid_,
      std::vector<Observation>& observations_);
    TransmissionModel (
      GlobalConfig* config_, 
      SpectralGrid* spectral_grid_,
      const size_t nb_grid_points_,
      const std::vector<std::string>& opacity_species_symbol,
      const std::vector<std::string>& opacity_species_folder);
    
    virtual ~TransmissionModel();

    virtual size_t parametersNumber() {
      return nb_total_param();};

    virtual bool calcModelCPU(
      const std::vector<double>& parameters, 
      std::vector<double>& spectrum, 
      std::vector<std::vector<double>>& spectrum_obs);
    
    virtual bool calcModelGPU(
      const std::vector<double>& parameter, 
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
      bool& delete_unused_files);

    virtual AtmosphereOutput getAtmosphereStructure(
      const std::vector<double>& physical_parameter,
      const std::vector<std::string>& species_symbols);
    
    virtual bool testModel(
      const std::vector<double>& parameters);

    std::vector<double> calcSpectrum(
      const double surface_gravity,
      const double planet_radius,
      const double radius_ratio,
      const std::vector<double>& pressure,
      const std::vector<double>& temperature,
      const std::vector<std::string>& species_symbol,
      const std::vector<std::vector<double>>& mixing_ratios,
      const std::vector<std::vector<double>>& cloud_optical_depth,
      const double use_variable_gravity);

  protected:
    Atmosphere atmosphere;
    OpacityCalculation opacity_calc;

    Temperature* temperature_profile = nullptr;
    std::vector<Chemistry*> chemistry;
    std::vector<CloudModel*> cloud_models;
    std::vector<Module*> modules;

    size_t nb_general_param = 0;
    size_t nb_total_chemistry_param = 0;
    size_t nb_temperature_param = 0;
    size_t nb_total_cloud_param = 0;
    size_t nb_total_modules_param = 0;

    size_t nb_total_param() {
        return nb_general_param 
          + nb_total_chemistry_param 
          + nb_temperature_param 
          + nb_total_cloud_param
          + nb_total_modules_param
          + nb_spectrum_modifier_param;};
    
    size_t nb_grid_points = 0;

    std::vector<std::vector<double>> cloud_extinction;
    float* cloud_extinction_gpu = nullptr;

    bool fit_mean_molecular_weight = false;
    bool fit_scale_height = false;
    bool use_variable_gravity = false;
    
    void initModules(
      const TransmissionModelConfig& model_config);
    
    std::vector<double> model_parameters;
    std::vector<double> chemistry_parameters;
    std::vector<double> cloud_parameters;
    std::vector<double> temperature_parameters;
    std::vector<double> module_parameters;
    std::vector<double> spectrum_modifier_parameters;
    
    void extractParameters(
      const std::vector<double>& parameters);
    
    bool calcAtmosphereStructure(
      const std::vector<double>& parameter);
    void setCloudProperties(
      const std::vector<std::vector<double>>& cloud_optical_depth);

    void postProcess(
      const TransmissionPostProcessConfig& post_process_config,
      const std::vector< std::vector<double> >& model_parameter,
      const size_t best_fit_model);

    void postProcessModel(
      const std::vector<double>& parameters, 
      std::vector<double>& temperature_profile, 
      std::vector<std::vector<double>>& mixing_ratios);

    void savePostProcessChemistry(
      const std::vector<std::vector<std::vector<double>>>& mixing_ratios, 
      const unsigned int species);
    void savePostProcessTemperatures(
      const std::vector<std::vector<double>>& temperature_profiles);
    
    void calcTransitDepthGPU(
      double* transit_radius_dev, 
      float* absorption_coeff_dev, 
      float* scattering_coeff_dev, 
      float* cloud_extinction_coeff_dev, 
      const Atmosphere& atmosphere, 
      const size_t nb_spectral_points, 
      const double radius_planet, 
      const double radius_star);

    void calcTransmissionSpectrum(
      const double bottom_radius, 
      const double star_radius, 
      std::vector<double>& spectrum);
    double calcTransitRadius(
      const unsigned int wavelength, 
      const double bottom_radius);
    double integrateEffectiveTangentHeight(
      const unsigned int wavelength, 
      const double bottom_radius);
    std::vector<double> tangentPathsTransmission(
      const unsigned int wavelength, 
      const double bottom_radius);
    double tangentOpticalDepth(
      const unsigned int tangent_altitude, 
      const unsigned int wavelength, 
      const double bottom_radius);
    double distanceToTangentCenter(
      const unsigned int tangent_altitude, 
      const unsigned int altitude, 
      const double bottom_radius);
};


}


#endif

