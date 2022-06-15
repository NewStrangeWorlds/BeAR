/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2020 Daniel Kitzmann
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


#ifndef _emission_h
#define _emission_h

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

#include "../forward_model.h"

#include "../atmosphere/atmosphere.h"
#include "../../cloud_model/cloud_model.h"
#include "../../chemistry/chemistry.h"
#include "../../temperature/temperature.h"
#include "../../transport_coeff/transport_coeff.h"
#include "../../radiative_transfer/radiative_transfer.h"


namespace helios {


//forward declaration
class Retrieval;


//this struct handles the Brown Dwarf config
//it will read in the corresponding parameter file
//and will then be used to create a model object
struct EmissionModelConfig{
  size_t nb_grid_points = 0;

  double atmos_boundaries[2] {0, 0};
  double atmos_top_pressure = 0;
  double atmos_bottom_pressure = 0;

  bool use_cloud_layer = false;

  std::string temperature_profile_model;
  std::vector<std::string> temperature_profile_parameters;

  std::string radiative_transfer_model;
  std::vector<std::string> radiative_transfer_parameters;

  std::vector<std::string> chemistry_model;
  std::vector<std::vector<std::string>> chemistry_parameters;

  std::string cloud_model;
  std::vector<std::string> cloud_model_parameters;

  std::vector<std::string> opacity_species_symbol;
  std::vector<std::string> opacity_species_folder;

  EmissionModelConfig (const std::string& folder_path);
  void readConfigFile(const std::string& file_name);
  void readChemistryConfig(std::fstream& file);
  void readOpacityConfig(std::fstream& file);
};




class EmissionModel : public ForwardModel{
  public:
    EmissionModel (Retrieval* retrieval_ptr, const EmissionModelConfig model_config);
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
    
    virtual bool testModel(const std::vector<double>& parameter, double* model_spectrum_gpu);
  protected:
    Retrieval* retrieval;
    TransportCoefficients transport_coeff;
    
    Atmosphere atmosphere;
    RadiativeTransfer* radiative_transfer = nullptr;
    Temperature* temperature_profile = nullptr;
    std::vector<Chemistry*> chemistry;
    CloudModel* cloud_model = nullptr;
    
    size_t nb_general_param = 0;
    size_t nb_total_chemistry_param = 0;
    size_t nb_temperature_param = 0;
    size_t nb_cloud_param = 0;

    size_t nb_total_param() {return nb_general_param + nb_total_chemistry_param + nb_temperature_param + nb_cloud_param;}
    
    size_t nb_grid_points = 0;

    std::vector< std::vector<double> > absorption_coeff;
    std::vector< std::vector<double> > scattering_coeff;

    std::vector< std::vector<double> > cloud_optical_depths;
    std::vector< std::vector<double> > cloud_single_scattering;
    std::vector< std::vector<double> > cloud_asym_param;

    //pointer to the array that holds the pointers to the coefficients on the GPU
    double* absorption_coeff_gpu = nullptr;
    double* scattering_coeff_dev = nullptr;

    double* cloud_optical_depths_dev = nullptr;
    double* cloud_single_scattering_dev = nullptr;
    double* cloud_asym_param_dev = nullptr;

    virtual void setPriors();
    void readPriorConfigFile(const std::string& file_name, std::vector<std::string>& prior_type, 
                                                           std::vector<std::string>& prior_description, 
                                                           std::vector<std::vector<double>>& prior_parameter);
    void initModules(const EmissionModelConfig& model_config);
    void initDeviceMemory();

    bool calcAtmosphereStructure(const std::vector<double>& parameter);
    double radiusDistanceScaling(const std::vector<double>& parameter);

    void postProcessSpectrum(std::vector<double>& model_spectrum, std::vector<double>& model_spectrum_bands);
    void postProcessSpectrumGPU(double* model_spectrum, double* model_spectrum_bands);

    void postProcessModel(const std::vector<double>& parameter, const std::vector<double>& model_spectrum_bands, 
                          std::vector<double>& temperature_profile, double& effective_temperature,
                          std::vector<std::vector<double>>& mixing_ratios);
    double postProcessEffectiveTemperature(const std::vector<double>& model_spectrum_bands, const double radius_distance_scaling);
    void savePostProcessChemistry(const std::vector<std::vector<std::vector<double>>>& mixing_ratios, const unsigned int species);
    void savePostProcessTemperatures(const std::vector<std::vector<double>>& temperature_profiles);
    void savePostProcessEffectiveTemperatures(const std::vector<double>& effective_temperatures);

    bool testCPUvsGPU(const std::vector<double>& parameter, double* model_spectrum_gpu);
};


}


#endif

