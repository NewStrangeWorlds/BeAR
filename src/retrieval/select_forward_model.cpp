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


#include "retrieval.h"

#include "../forward_model/forward_model.h"
#include "../additional/exceptions.h"

//the header files for all forward models
#include "../forward_model/emission/emission.h"
#include "../forward_model/secondary_eclipse/secondary_eclipse.h"
#include "../forward_model/transmission/transmission.h"
#include "../forward_model/flat_line/flat_line.h"
#include "../forward_model/secondary_eclipse_bb/secondary_eclipse_bb.h"


namespace bear{

//Selects and initialises the forward model based on the option found in retrieval.config
//Exits with an error if the selected forward model is unkown
ForwardModel* Retrieval::selectForwardModel(
  const std::string model_description,
  GenericConfig* model_config)
{
  if (model_description == "emission" || model_description == "Emission" || model_description == "em")
  {
    EmissionModel* model = new EmissionModel(
      EmissionModelConfig (config->retrieval_folder_path),
      config,
      &spectral_grid,
      observations);

    return model;
  }


  if (model_description == "secondary_eclipse" || model_description == "Secondary_eclipse" || model_description == "se")
  {
    SecondaryEclipseModel* model = new SecondaryEclipseModel(
      SecondaryEclipseConfig (config->retrieval_folder_path),
      config,
      &spectral_grid,
      observations);

    return model;
  }


  if (model_description == "transmission" || model_description == "Transmission" || model_description == "trans")
  {
    if (model_config == nullptr)
    {
      TransmissionModel* model = new TransmissionModel(
        TransmissionModelConfig (config->retrieval_folder_path),
        config,
        &spectral_grid,
        observations);

      return model;
    }
    else
    {
      TransmissionModelConfig* c = dynamic_cast<TransmissionModelConfig*>(model_config);

      TransmissionModel* model = new TransmissionModel(
        *c,
        config,
        &spectral_grid,
        observations);

      return model;
    }
  }


  if (model_description == "flat_line" || model_description == "Flat_line" || model_description == "fl")
  {
    FlatLine* model = new FlatLine(
      config,
      &spectral_grid,
      observations);

    return model;
  }


  if (model_description == "secondary_eclipse_bb" || model_description == "Secondary_eclipse_bb" || model_description == "se_bb")
  {
    SecondaryEclipseBlackBodyModel* model = new SecondaryEclipseBlackBodyModel(
      SecondaryEclipseBlackBodyConfig (config->retrieval_folder_path),
      config,
      &spectral_grid,
      observations);

    return model;
  }


  std::string error_message = "Unkown forward model found in retrieval config file: " + model_description + "\n";
  throw InvalidInput(std::string ("retrieval.config"), error_message);

  return nullptr;
}


}
