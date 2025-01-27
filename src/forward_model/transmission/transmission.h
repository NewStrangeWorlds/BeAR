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
#include "../../retrieval/priors.h"
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
    double atmos_top_pressure = 0;
    double atmos_bottom_pressure = 0;

    bool fit_mean_molecular_weight = false;
    bool fit_scale_height = false;
    bool use_variable_gravity = false;

    bool use_cloud_model = false;
    bool use_optional_modules = false;

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
    virtual void readConfigFile(const std::string& file_name);
};


class TransmissionPostProcessConfig : public GenericConfig{
  public:
    std::vector<chemical_species_id> species_to_save;

    bool save_temperatures = false;
    bool save_spectra = true;

    bool delete_sampler_files = false;

    TransmissionPostProcessConfig (const std::string& folder_path);
    void readConfigFile(const std::string& file_name);
};


class TransmissionModel : public ForwardModel{
  public:
    TransmissionModel (
      const TransmissionModelConfig model_config,
      Priors* priors_,
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
    
    virtual bool calcModel(
      const std::vector<double>& parameter, 
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
          + nb_spectrum_modifier_param;
      }
    
    size_t nb_grid_points = 0;

    std::vector<std::vector<double>> cloud_extinction;
    double* cloud_extinction_gpu = nullptr;

    bool fit_mean_molecular_weight = false;
    bool fit_scale_height = false;
    bool use_variable_gravity = false;
    
    virtual void setPriors(Priors* priors);
    
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
      double* absorption_coeff_dev, 
      double* scattering_coeff_dev, 
      double* cloud_extinction_coeff_dev, 
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

