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


#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "secondary_eclipse.h"

#include "../../additional/exceptions.h"
#include "../../chemistry/select_chemistry.h"
#include "../stellar_spectrum/select_stellar_model.h"
#include "../../radiative_transfer/select_radiative_transfer.h"
#include "../../temperature/select_temperature_profile.h"
#include "../../cloud_model/select_cloud_model.h"
#include "../../CUDA_kernels/data_management_kernels.h"


namespace bear{


//initialises the varous modules of the forward model
void SecondaryEclipseModel::initModules(const SecondaryEclipseConfig& model_config)
{
  radiative_transfer = selectRadiativeTransfer(
    model_config.radiative_transfer_model, 
    model_config.radiative_transfer_parameters, 
    model_config.nb_grid_points, 
    config, 
    spectral_grid);


  stellar_model = selectStellarModel(
    model_config.stellar_spectrum_model,
    model_config.stellar_model_parameters,
    spectral_grid);

  nb_stellar_param = stellar_model->nbParameters();


  chemistry.assign(model_config.chemistry_model.size(), nullptr);

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
    CloudModel* model = selectCloudModel(
      model_config.cloud_model[i], 
      model_config.cloud_model_parameters[i]);
    
    if (model != nullptr)
      cloud_models.push_back(model);
  }
  
  //count the total number of free parameters for the cloud modules
  nb_total_cloud_param = 0;
  
  for (auto & i : cloud_models)
    nb_total_cloud_param += i->nbParameters();
}



}

