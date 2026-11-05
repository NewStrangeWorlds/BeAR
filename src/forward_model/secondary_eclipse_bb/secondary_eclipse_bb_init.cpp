// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "secondary_eclipse_bb.h"


#include "../../additional/exceptions.h"
#include "../../CUDA_kernels/data_management_kernels.h"
#include "../stellar_spectrum/select_stellar_model.h"


namespace bear{


//initialises the varous modules of the forward model
void OccultationBlackBodyModel::initModules(
  const OccultationBlackBodyConfig& model_config)
{
  stellar_model = selectStellarModel(
    model_config.stellar_spectrum_model,
    model_config.stellar_model_parameters,
    spectral_grid);

  nb_stellar_param = stellar_model->nbParameters();


  //Assemble the ordered list of parameter names. The order here MUST match the
  //slicing in OccultationBlackBodyModel::extractParameters:
  //  general | stellar | spectrum modifier
  parameter_names.clear();

  parameter_names.push_back("planet_temperature");
  parameter_names.push_back("radius_ratio");

  for (auto & name : stellar_model->parameterNames())
    parameter_names.push_back(name);

  for (size_t i=0; i<nb_spectrum_modifier_param; ++i)
    parameter_names.push_back("spectrum_shift_" + std::to_string(i+1));
}



}

