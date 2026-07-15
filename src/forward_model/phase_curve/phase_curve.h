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


#ifndef _phase_curve_h
#define _phase_curve_h

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <memory>

#include "../forward_model.h"
#include "../generic_config.h"
#include "../modules/module.h"

#include "../../config/global_config.h"
#include "../../spectral_grid/spectral_grid.h"
#include "../../observations/observations.h"
#include "../atmosphere/atmosphere.h"
#include "../../chemistry/chemistry.h"
#include "../../cloud_model/cloud_model.h"
#include "../../temperature/temperature.h"
#include "../../transport_coeff/transport_coeff.h"
#include "../../transport_coeff/opacity_calc.h"
#include "../../radiative_transfer/radiative_transfer.h"
#include "../stellar_spectrum/stellar_spectrum.h"


namespace bear {


struct PhaseCurveConfig : public GenericConfig{
  size_t nb_grid_points = 0;

  std::vector<double> atmos_boundaries = {0, 0};

  std::string temperature_profile_model;
  std::vector<std::string> temperature_profile_parameters;

  std::string stellar_spectrum_model;
  std::vector<std::string> stellar_model_parameters;

  std::string radiative_transfer_model;
  std::vector<std::string> radiative_transfer_parameters;

  std::vector<std::string> chemistry_model;
  std::vector<std::vector<std::string>> chemistry_parameters;

  std::vector<std::string> cloud_model;
  std::vector<std::vector<std::string>> cloud_model_parameters;

  std::vector<std::string> modules;
  std::vector<std::vector<std::string>> modules_parameters;

  std::vector<std::string> opacity_species_symbol;
  std::vector<std::string> opacity_species_folder;

  std::vector<std::string> opacity_species_symbol_highres;
  std::vector<std::string> opacity_species_folder_highres;

  double highres_stellar_smooth_sigma = 0.0;

  PhaseCurveConfig (
    const std::string& folder_path,
    const std::string& file_name = "forward_model.toml");
  PhaseCurveConfig (
    const int nb_grid_points_,
    const double atmos_bottom_pressure_,
    const double atmos_top_pressure_,
    const std::string temperature_profile_model_,
    const std::vector<std::string>& temperature_profile_parameters_,
    const std::string radiative_transfer_model_,
    const std::vector<std::string>& radiative_transfer_parameters_,
    const std::vector<std::string>& chemistry_model_,
    const std::vector<std::vector<std::string>>& chemistry_parameters_,
    const std::vector<std::string>& opacity_species_symbol_,
    const std::vector<std::string>& opacity_species_folder_);
  PhaseCurveConfig (
    const int nb_grid_points_,
    const double atmos_bottom_pressure_,
    const double atmos_top_pressure_,
    const std::string temperature_profile_model_,
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



class PhaseCurvePostProcessConfig : public GenericConfig{
  public:
    bool save_temperatures = true;
    bool save_spectra = true;
    bool save_contribution_functions = false;
    bool delete_sampler_files = false;
    std::vector<chemical_species_id> species_to_save;

    PhaseCurvePostProcessConfig (
      const std::string& folder_path,
      const std::string& file_name = "post_process.toml");
    PhaseCurvePostProcessConfig (
      const bool save_temperatures_,
      const bool save_spectra_,
      const bool save_contribution_functions_,
      const std::vector<std::string>& species_to_save_);

    void readConfigFile(const std::string& file_name);
};




class PhaseCurveModel : public ForwardModel{
  public:
    PhaseCurveModel (
      const PhaseCurveConfig model_config,
      GlobalConfig* config_,
      SpectralGrid* spectral_grid_,
      std::vector<Observation>& observations_);

    virtual ~PhaseCurveModel();

    virtual size_t parametersNumber() {
      //parameter_names is the source of truth once assembled in initModules;
      //the lightweight test constructor leaves it empty and falls back.
      return parameter_names.empty() ? nb_total_param() : parameter_names.size();};

    virtual bool calcModelCPU(
      const std::vector<double>& parameter,
      std::vector<double>& spectrum,
      std::vector<std::vector<double>>& spectrum_obs);
    virtual bool calcModelGPU(
      const std::vector<double>& parameters,
      float* spectrum,
      std::vector<float*>& spectrum_obs);

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

    std::vector<double> calcSpectrum(
      const double surface_gravity,
      const std::vector<double>& pressure,
      const std::vector<double>& temperature,
      const std::vector<std::string>& species_symbol,
      const std::vector<std::vector<double>>& mixing_ratios,
      const std::vector<std::vector<double>>& cloud_optical_depth);
    virtual void setHighResGrid(SpectralGrid* grid) override;

  protected:
    Atmosphere atmosphere;
    OpacityCalculation opacity_calc;
    std::unique_ptr<OpacityCalculation> opacity_calc_highres;

    std::unique_ptr<RadiativeTransfer> radiative_transfer;
    std::unique_ptr<RadiativeTransfer> radiative_transfer_highres;
    std::unique_ptr<Temperature> temperature_profile;
    std::unique_ptr<StellarSpectrumModel> stellar_model;
    std::unique_ptr<StellarSpectrumModel> stellar_model_highres_;
    std::vector<std::unique_ptr<Chemistry>> chemistry;
    std::vector<std::unique_ptr<CloudModel>> cloud_models;
    std::vector<std::unique_ptr<Module>> modules;

    std::vector<std::string> opacity_species_symbol_;
    std::vector<std::string> opacity_species_folder_;
    std::vector<std::string> opacity_species_symbol_highres_;
    std::vector<std::string> opacity_species_folder_highres_;
    std::string radiative_transfer_model_;
    std::vector<std::string> radiative_transfer_parameters_;
    std::string stellar_spectrum_model_name_;
    std::vector<std::string> stellar_model_parameters_names_;
    double highres_stellar_smooth_sigma_ = 0.0;

    std::vector<size_t> modules_lowres_idx;
    std::vector<size_t> modules_highres_idx;

    size_t nb_grid_points = 0;
    size_t nb_general_param = 0;
    size_t nb_stellar_param = 0;
    size_t nb_total_chemistry_param = 0;
    size_t nb_temperature_param = 0;
    size_t nb_total_cloud_param = 0;
    size_t nb_total_modules_param = 0;

    size_t nb_total_param() {
      return nb_general_param
             + nb_stellar_param
             + nb_total_chemistry_param
             + nb_temperature_param
             + nb_total_cloud_param
             + nb_total_modules_param
             + nb_spectrum_modifier_param;
    }

    void initModules(const PhaseCurveConfig& model_config);

    float* stellar_flux_highres_gpu_ = nullptr;

    // Cached stellar spectrum for per-pixel Doppler correction in high-res kernels.
    std::vector<double> stellar_spectrum_cpu_;

    const std::vector<double>& stellarSpectrumCPU() const override
      { return stellar_spectrum_cpu_; }
    const float* stellarSpectrumGPU() const override
      { return stellar_flux_highres_gpu_; }

    void normaliseFpFsGPU(
      float*       planet_spectrum,
      const float* stellar_spectrum,
      const int    nb_points,
      const float  radius_ratio_squared);

    void calcOccultationLowResGPU(
      float*       output,
      const float* planet_spectrum,
      const float* stellar_spectrum,
      const int    nb_points,
      const float  radius_ratio_squared);

    std::vector<double> model_parameters;
    std::vector<double> stellar_parameters;
    std::vector<double> chemistry_parameters;
    std::vector<double> cloud_parameters;
    std::vector<double> temperature_parameters;
    std::vector<double> module_parameters;
    std::vector<double> spectrum_modifier_parameters;

    void extractParameters(
      const std::vector<double>& parameters);

    bool calcAtmosphereStructure(const std::vector<double>& parameter);

    void setCloudProperties(const std::vector<std::vector<double>>& cloud_optical_depth);

    void postProcess(
      const PhaseCurvePostProcessConfig& post_process_config,
      const std::vector< std::vector<double> >& model_parameter,
      const size_t best_fit_model);
    void postProcessModel(
      const std::vector<double>& parameter,
      std::vector<double>& temperature_profile,
      std::vector<std::vector<double>>& mixing_ratios);
    void postProcessContributionFunctions(
      const std::vector<double>& model_parameter);
    void saveContributionFunctions(
      std::vector< std::vector<double>>& contribution_function,
      const size_t observation_index);
    void savePostProcessChemistry(
      const std::vector<std::vector<std::vector<double>>>& mixing_ratios,
      const unsigned int species);
    void savePostProcessTemperatures(
      const std::vector<std::vector<double>>& temperature_profiles);
};


}


#endif


