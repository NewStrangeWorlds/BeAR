/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2022 Daniel Kitzmann
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


#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "secondary_eclipse_bb.h"


#include "../../additional/exceptions.h"
#include "../../CUDA_kernels/data_management_kernels.h"
#include "../stellar_spectrum/select_stellar_model.h"


namespace helios{


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

