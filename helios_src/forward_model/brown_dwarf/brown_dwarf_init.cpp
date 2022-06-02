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


#include "brown_dwarf.h"


#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>


#include "../../additional/exceptions.h"
#include "../../retrieval/retrieval.h"

#include "../../chemistry/select_chemistry.h"
#include "../../temperature/select_temperature_profile.h"
#include "../../radiative_transfer/select_radiative_transfer.h"


namespace helios{


//initialise radiative transfer model
void BrownDwarfModel::initRadiativeTransfer(const BrownDwarfConfig& model_config)
{

  radiative_transfer = selectRadiativeTransfer(model_config.radiative_transfer_model, 
                                               model_config.radiative_transfer_parameters, 
                                               model_config.nb_grid_points, 
                                               retrieval->config, &retrieval->spectral_grid);

}



//select and initialise the chemistry models
void BrownDwarfModel::initChemistry(const BrownDwarfConfig& model_config)
{
  chemistry.assign(model_config.chemistry_model.size(), nullptr);

  for (size_t i=0; i<model_config.chemistry_model.size(); ++i)
    chemistry[i] = selectChemistryModule(model_config.chemistry_model[i], model_config.chemistry_parameters[i], retrieval->config, model_config.atmos_boundaries);

  nb_total_chemistry_param = 0;

  for (auto & i : chemistry)
    nb_total_chemistry_param += i->nbParameters(); 
}



//select and initialise the chemistry models
void BrownDwarfModel::initTemperature(const BrownDwarfConfig& model_config)
{
   temperature_profile = selectTemperatureProfile(model_config.temperature_profile_model, 
                                                  model_config.temperature_profile_parameters, 
                                                  model_config.atmos_boundaries);

}



}

