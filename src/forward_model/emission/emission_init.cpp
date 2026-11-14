// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "emission.h"

#include "../../additional/exceptions.h"
#include "../../chemistry/select_chemistry.h"
#include "../../temperature/select_temperature_profile.h"
#include "../../radiative_transfer/select_radiative_transfer.h"
#include "../../cloud_model/select_cloud_model.h"


namespace bear{


//initialises the varous modules of the forward model
void EmissionModel::initModules(const EmissionModelConfig& model_config)
{
  radiative_transfer = selectRadiativeTransfer(
    model_config.radiative_transfer_model, 
    model_config.radiative_transfer_parameters, 
    model_config.nb_grid_points, 
    config, 
    spectral_grid);


  chemistry.resize(model_config.chemistry_model.size());

  for (size_t i=0; i<model_config.chemistry_model.size(); ++i)
    chemistry[i] = selectChemistryModule(
      model_config.chemistry_model[i],
      model_config.chemistry_parameters[i],
      config,
      model_config.atmos_boundaries);

  //count the total number of free parameters for the chemistry modules
  nb_total_chemistry_param = 0;

  for (auto & i : chemistry)
    nb_total_chemistry_param += i->nbParameters(); 


  temperature_profile = selectTemperatureProfile(
    model_config.temperature_profile_model, 
    model_config.temperature_profile_parameters, 
    model_config.atmos_boundaries);

  nb_temperature_param = temperature_profile->nbParameters();


  for (size_t i=0; i<model_config.cloud_model.size(); ++i)
  {
    auto model = selectCloudModel(
      model_config.cloud_model[i],
      model_config.cloud_model_parameters[i],
      config);

    if (model != nullptr)
      cloud_models.push_back(std::move(model));
  }
  
  //count the total number of free parameters for the cloud modules
  nb_total_cloud_param = 0;

  for (auto & i : cloud_models)
    nb_total_cloud_param += i->nbParameters();


  //Assemble the ordered list of parameter names for this model's parameter
  //block. The order here MUST match the slicing in EmissionModel::extractParameters:
  //  general | chemistry | temperature | cloud | spectrum modifier
  parameter_names.clear();

  parameter_names.push_back("log_g");
  parameter_names.push_back("scaling_factor");
  parameter_names.push_back("distance");

  for (auto & i : chemistry)
    for (auto & name : i->parameterNames())
      parameter_names.push_back(name);

  for (auto & name : temperature_profile->parameterNames())
    parameter_names.push_back(name);

  for (size_t c=0; c<cloud_models.size(); ++c)
    appendParameterNames(cloud_models[c]->parameterNames(), c, cloud_models.size());

  for (size_t i=0; i<nb_spectrum_modifier_param; ++i)
    parameter_names.push_back("spectrum_shift_" + std::to_string(i+1));
}


}

