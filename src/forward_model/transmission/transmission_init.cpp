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

#include "transmission.h"

#include "../../chemistry/select_chemistry.h"
#include "../../temperature/select_temperature_profile.h"
#include "../../cloud_model/select_cloud_model.h"
#include "../modules/select_module.h"
#include "../modules/velocity_broadening/velocity_broadening.h"
#include "../../transport_coeff/opacity_calc.h"


namespace bear{


//initialises the varous modules of the forward model
void TransmissionModel::initModules(const TransmissionModelConfig& model_config)
{ 
  chemistry.resize(model_config.chemistry_model.size());

  for (size_t i=0; i<model_config.chemistry_model.size(); ++i)
    chemistry[i] = selectChemistryModule(
      model_config.chemistry_model[i],
      model_config.chemistry_parameters[i],
      config,
      model_config.atmos_boundaries);

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
      model_config.cloud_model_parameters[i]);

    if (model != nullptr)
      cloud_models.push_back(std::move(model));
  }
  
  nb_total_cloud_param = 0;
  
  for (auto & i : cloud_models)
    nb_total_cloud_param += i->nbParameters();

  if (model_config.modules.size() > 0)
  {
    for (size_t i=0; i<model_config.modules.size(); ++i)
    {
      auto module = selectModule(
        model_config.modules[i],
        model_config.modules_parameters[i],
        spectral_grid);

      if (module != nullptr)
        modules.push_back(std::move(module));
    }
  
    nb_total_modules_param = 0;
  
    for (auto & i : modules)
      nb_total_modules_param += i->nbParameters();
    }
}




void TransmissionModel::setHighResGrid(SpectralGrid* grid)
{
  ForwardModel::setHighResGrid(grid);

  if (!spectral_grid_highres) return;

  std::cout << "Initialising high-res spectral grid for transmission model\n";

  // Create a second OpacityCalculation on the high-res grid,
  // sharing the same Atmosphere object
  opacity_calc_highres = std::make_unique<OpacityCalculation>(
    config,
    spectral_grid_highres,
    &atmosphere,
    opacity_species_symbol_,
    opacity_species_folder_,
    config->use_gpu,
    false);  // no clouds on high-res grid

  // Classify modules: velocity broadening -> highres, everything else -> lowres
  for (size_t i = 0; i < modules.size(); ++i)
  {
    auto* vb = dynamic_cast<VelocityBroadening*>(modules[i].get());

    if (vb != nullptr)
    {
      modules_highres_idx.push_back(i);
      vb->setSpectralGrid(spectral_grid_highres);
    }
    else
      modules_lowres_idx.push_back(i);
  }

  // Allocate high-res spectrum GPU buffer
  if (config->use_gpu)
    allocateOnDevice(
      spectrum_highres_gpu_,
      spectral_grid_highres->nbSpectralPoints());

  std::cout << "High-res grid: " << spectral_grid_highres->nbSpectralPoints()
            << " spectral points\n";
  std::cout << "Modules on low-res path: " << modules_lowres_idx.size()
            << ", on high-res path: " << modules_highres_idx.size() << "\n";
}


}