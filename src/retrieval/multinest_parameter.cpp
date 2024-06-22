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


#include "multinest_parameter.h"

#include "../config/global_config.h"

#include <string>
#include <iostream>
#include <fstream>
#include <cstring>



namespace bear{


//intialise the parameters from the global config
MultinestParameter::MultinestParameter(GlobalConfig* config)
{
  is = config->multinest_ins;
  mmodal = config->multinest_mode_sep;
  nlive = config->multinest_nb_living_points;
  efr = config->multinest_efficiency;
  maxiter = config->multinest_nb_iterations;
  resume = config->multinest_resume;
  fb = config->multinest_feedback;
  std::strcpy(root, config->retrieval_folder_path.c_str());
}



//not implemented yet :(
void MultinestParameter::loadParameterFile(const std::string& file_name)
{


}



}
