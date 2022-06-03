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


#include "retrieval.h"
#include "../forward_model/forward_model.h"
#include "../additional/exceptions.h"


//the header files for all forward models
#include "../forward_model/brown_dwarf/brown_dwarf.h"
#include "../forward_model/secondary_eclipse/secondary_eclipse.h"



namespace helios{


//Selects and initialises the forward model based on the option found in retrieval.config
//Exits with an error if the selected forward model is unkown
ForwardModel* Retrieval::selectForwardModel(const std::string model_description)
{
  if (model_description == "emission" || model_description == "Emission" || model_description == "em")
  {
    BrownDwarfModel* model = new BrownDwarfModel(this, BrownDwarfConfig (config->retrieval_folder_path));

    return model;
  }


  if (model_description == "secondary_eclipse" || model_description == "Secondary_eclipse" || model_description == "se")
  {
    SecondaryEclipseModel* model = new SecondaryEclipseModel(this, SecondaryEclipseConfig (config->retrieval_folder_path));

    return model;
  }


  std::string error_message = "Unkown forward model found in retrieval config file: " + model_description + "\n";
  throw ExceptionInvalidInput(std::string ("retrieval.config"), error_message);

  return nullptr;
}


}
