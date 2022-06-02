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


#ifndef _brown_dwarf_h
#define _brown_dwarf_h

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

#include "../forward_model.h"

#include "../atmosphere/atmosphere.h"

#include "../../chemistry/chemistry.h"
#include "../../temperature/temperature.h"
#include "../../transport_coeff/transport_coeff.h"
#include "../../radiative_transfer/discrete_ordinate.h"
#include "../../radiative_transfer/short_characteristics.h"


#include "../../radiative_transfer/radiative_transfer.h"


namespace helios {


//forward declaration
class Retrieval;


//this struct handles the Brown Dwarf config
//it will read in the corresponding parameter file
//and will then be used to create a model object
struct BrownDwarfConfig{
  size_t nb_grid_points = 0;

  double atmos_boundaries[2] {0, 0};
  double atmos_top_pressure = 0;
  double atmos_bottom_pressure = 0;
  
  size_t nb_temperature_elements = 0;
  size_t temperature_poly_degree = 0;

  bool use_cloud_layer = false;
  size_t radiative_transfer_model = 0;
  
  std::vector<std::string> chemistry_model;
  std::vector<std::vector<std::string>> chemistry_parameters;
  
  std::vector<std::string> opacity_species_symbol;
  std::vector<std::string> opacity_species_folder;
  
  BrownDwarfConfig (const std::string& folder_path);
  void readConfigFile(const std::string& file_name);
  void readChemistryConfig(std::fstream& file);
  void readOpacityConfig(std::fstream& file);
};




class BrownDwarfModel : public ForwardModel{
  public:
    BrownDwarfModel (Retrieval* retrieval_ptr, const BrownDwarfConfig model_config);
    virtual ~BrownDwarfModel();
    virtual bool calcModel(const std::vector<double>& parameter, std::vector<double>& spectrum);
    virtual bool calcModelGPU(const std::vector<double>& parameter, double* model_spectrum, double* model_spectrum_bands);
    
    virtual void postProcess(const std::vector< std::vector<double> >& model_parameter, 
                             const std::vector< std::vector<double> >& model_spectrum_bands,
                             const size_t best_fit_model);
  protected:
    Retrieval* retrieval;
    TransportCoefficients transport_coeff;
    RadiativeTransfer* radiative_transfer = nullptr;

    Atmosphere atmosphere;
    Temperature* temperature_profile = nullptr;
    std::vector<Chemistry*> chemistry;


    size_t nb_grid_points = 0;
    size_t nb_general_param = 0;
    size_t nb_total_chemistry_param = 0;
    size_t nb_cloud_param = 0;

    size_t nb_total_param() {return nb_general_param + nb_total_chemistry_param + temperature_profile->nbParameters() + nb_cloud_param;}

    double radius_distance_scaling = 0;

    std::vector< std::vector<double> > absorption_coeff;
    std::vector< std::vector<double> > scattering_coeff;

    bool use_cloud_layer = false;
    std::vector<double> cloud_optical_depths;

    //pointer to the array that holds the pointers to the coefficients on the GPU
    double* absorption_coeff_gpu = nullptr;
    
    virtual void setPriors();
    void readPriorConfigFile(const std::string& file_name, std::vector<std::string>& prior_type, 
                                                           std::vector<std::string>& prior_description, 
                                                           std::vector<std::vector<double>>& prior_parameter);
    void initChemistry(const BrownDwarfConfig& model_config);
    void initRadiativeTransfer(const BrownDwarfConfig& model_config);
    void initTemperature(const BrownDwarfConfig& model_config);

    bool calcAtmosphereStructure(const std::vector<double>& parameter);
    double radiusDistanceScaling(const double distance, const double radius, const double scaling_f);

    void calcGreyCloudLayer(const std::vector<double>& cloud_parameters);
    void calcCloudPosition(const double top_pressure, const double bottom_pressure, unsigned int& top_index, unsigned int& bottom_index);

    void postProcessModel(const std::vector<double>& parameter, const std::vector<double>& model_spectrum_bands, 
                          std::vector<double>& temperature_profile, double& effective_temperature,
                          std::vector<std::vector<double>>& mixing_ratios);
    double postProcessEffectiveTemperature(const std::vector<double>& model_spectrum_bands);
    void savePostProcessChemistry(const std::vector<std::vector<std::vector<double>>>& mixing_ratios, const unsigned int species);
    void savePostProcessTemperatures(const std::vector<std::vector<double>>& temperature_profiles);
    void savePostProcessEffectiveTemperatures(const std::vector<double>& effective_temperatures);
};


}


#endif

