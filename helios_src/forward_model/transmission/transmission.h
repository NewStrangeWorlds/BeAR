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


#ifndef _transmission_h
#define _transmission_h


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
#include "../../transport_coeff/opacity_calc.h"


namespace helios {


constexpr double transmission_optical_depth_cutoff = 10;
constexpr double transmission_cutoff = 4.5399929e-5; //exp(-tau)


//this struct handles the Transmission Spectrum config
//it will read in the corresponding parameter file
//and will then be used to create a model object
struct TransmissionModelConfig{
  size_t nb_grid_points = 0;

  double atmos_boundaries[2] {0, 0};
  double atmos_top_pressure = 0;
  double atmos_bottom_pressure = 0;

  bool use_cloud_model = false;

  std::string temperature_profile_model;
  std::vector<std::string> temperature_profile_parameters;

  std::vector<std::string> chemistry_model;
  std::vector<std::vector<std::string>> chemistry_parameters;

  std::vector<std::string> cloud_model;
  std::vector<std::vector<std::string>> cloud_model_parameters;

  std::vector<std::string> opacity_species_symbol;
  std::vector<std::string> opacity_species_folder;

  TransmissionModelConfig (const std::string& folder_path);
  void readConfigFile(const std::string& file_name);
  void readCloudConfig(std::fstream& file);
  void readChemistryConfig(std::fstream& file);
  void readOpacityConfig(std::fstream& file);
};



class TransmissionModel : public ForwardModel{
  public:
    TransmissionModel (
      const TransmissionModelConfig model_config,
      Priors* priors_,
      GlobalConfig* config_, 
      SpectralGrid* spectral_grid_,
      std::vector<Observation>& observations_);
    virtual ~TransmissionModel();
    
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

    Atmosphere atmosphere;
    OpacityCalculation opacity_calc;

    Temperature* temperature_profile = nullptr;
    std::vector<Chemistry*> chemistry;
    CloudModel* cloud_model = nullptr;
    std::vector<CloudModel*> cloud_models;

    std::vector<Observation>& observations;
    size_t nb_observation_points = 0;

    size_t nb_general_param = 0;
    size_t nb_total_chemistry_param = 0;
    size_t nb_temperature_param = 0;
    size_t nb_total_cloud_param = 0;

    size_t nb_total_param() {
        return nb_general_param + nb_total_chemistry_param + nb_temperature_param + nb_total_cloud_param;
      }
    
    size_t nb_grid_points = 0;

    std::vector<std::vector<double>> cloud_extinction;
    double* cloud_extinction_gpu = nullptr;

    virtual void setPriors(Priors* priors);

    void readPriorConfigFile(
      const std::string& file_name, 
      std::vector<std::string>& prior_type, 
      std::vector<std::string>& prior_description, 
      std::vector<std::vector<double>>& prior_parameter);
    
    void initModules(const TransmissionModelConfig& model_config);

    bool calcAtmosphereStructure(const std::vector<double>& parameter);

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

    bool testCPUvsGPU(
      const std::vector<double>& parameter, 
      double* model_spectrum_gpu);
};


}


#endif

