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
#include "../forward_model/phase_curve/phase_curve.h"

#include <memory>


namespace bear{

//Selects and initialises the forward model based on the option found in retrieval.config
//Exits with an error if the selected forward model is unkown
std::unique_ptr<ForwardModel> Retrieval::selectForwardModel(
  const std::string model_description,
  GenericConfig* model_config)
{
  if (model_description == "emission" || model_description == "Emission" || model_description == "em")
  {
    if (model_config == nullptr)
    {
      return std::make_unique<EmissionModel>(
        EmissionModelConfig (config->retrieval_folder_path),
        config,
        &spectral_grid,
        observations);
    }
    else
    {
      EmissionModelConfig* c = dynamic_cast<EmissionModelConfig*>(model_config);

      return std::make_unique<EmissionModel>(
        *c,
        config,
        &spectral_grid,
        observations);
    }
  }


  if (model_description == "secondary_eclipse" || model_description == "Secondary_eclipse" || model_description == "se")
  {
    if (model_config == nullptr)
    {
      return std::make_unique<OccultationModel>(
        OccultationConfig (config->retrieval_folder_path),
        config,
        &spectral_grid,
        observations);
    }
    {
      OccultationConfig* c = dynamic_cast<OccultationConfig*>(model_config);

      return std::make_unique<OccultationModel>(
        *c,
        config,
        &spectral_grid,
        observations);
    }
  }


  if (model_description == "transmission" || model_description == "Transmission" || model_description == "trans")
  {
    if (model_config == nullptr)
    {
      return std::make_unique<TransmissionModel>(
        TransmissionModelConfig (config->retrieval_folder_path),
        config,
        &spectral_grid,
        observations);
    }
    else
    {
      TransmissionModelConfig* c = dynamic_cast<TransmissionModelConfig*>(model_config);

      return std::make_unique<TransmissionModel>(
        *c,
        config,
        &spectral_grid,
        observations);
    }
  }


  if (model_description == "flat_line" || model_description == "Flat_line" || model_description == "fl")
  {
    return std::make_unique<FlatLine>(
      config,
      &spectral_grid,
      observations);
  }


  if (model_description == "secondary_eclipse_bb" || model_description == "Secondary_eclipse_bb" || model_description == "se_bb")
  {
    if (model_config == nullptr)
    {
      return std::make_unique<OccultationBlackBodyModel>(
        OccultationBlackBodyConfig (config->retrieval_folder_path),
        config,
        &spectral_grid,
        observations);
    }
    else
    {
      OccultationBlackBodyConfig* c = dynamic_cast<OccultationBlackBodyConfig*>(model_config);

      return std::make_unique<OccultationBlackBodyModel>(
        *c,
        config,
        &spectral_grid,
        observations);
    }
  }


  if (model_description == "phase_curve" || model_description == "Phase_curve" || model_description == "pc")
  {
    if (model_config == nullptr)
    {
      return std::make_unique<PhaseCurveModel>(
        PhaseCurveConfig(config->retrieval_folder_path),
        config,
        &spectral_grid,
        observations);
    }
    else
    {
      PhaseCurveConfig* c = dynamic_cast<PhaseCurveConfig*>(model_config);

      return std::make_unique<PhaseCurveModel>(
        *c,
        config,
        &spectral_grid,
        observations);
    }
  }


  std::string error_message = "Unkown forward model found in retrieval config file: " + model_description + "\n";
  throw InvalidInput(std::string ("retrieval.config"), error_message);

  return nullptr;
}


}
