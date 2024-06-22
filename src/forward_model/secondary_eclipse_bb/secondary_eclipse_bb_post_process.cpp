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
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "secondary_eclipse_bb.h"

#include "../../chemistry/chem_species.h"
#include "../../CUDA_kernels/data_management_kernels.h"
#include "../../CUDA_kernels/contribution_function_kernels.h"
#include "../atmosphere/atmosphere.h"


namespace helios{


//calls the model specific posterior calculations
void SecondaryEclipseBlackBodyModel::postProcess(
  const std::vector< std::vector<double> >& model_parameter, 
  const std::vector< std::vector<double> >& model_spectrum_bands, 
  const size_t best_fit_model)
{
  //nothing to do here 

}


std::vector<double> SecondaryEclipseBlackBodyModel::convertSpectrumToModel(const std::vector<double>& spectrum)
{
  //the high-res spectrum is already a secondary eclipse depth in ppm
  std::vector<double> model_spectrum = spectrum;

  return model_spectrum;
}



}

