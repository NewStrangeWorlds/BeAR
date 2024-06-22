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

#include "secondary_eclipse_bb.h"


#include "../../additional/exceptions.h"
#include "../../CUDA_kernels/data_management_kernels.h"
#include "../stellar_spectrum/select_stellar_model.h"


namespace bear{


//initialises the varous modules of the forward model
void SecondaryEclipseBlackBodyModel::initModules(
  const SecondaryEclipseBlackBodyConfig& model_config)
{
  stellar_model = selectStellarModel(
    model_config.stellar_spectrum_model,
    model_config.stellar_model_parameters,
    spectral_grid);

  nb_stellar_param = stellar_model->nbParameters();
}



}

