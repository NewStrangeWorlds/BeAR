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


#ifndef _select_chemistry_h
#define _select_chemistry_h


#include "chemistry.h"

#include "fastchem_chemistry.h"
#include "isoprofile_chemistry.h"
#include "isoprofile_clr_chemistry.h"
#include "free_chemistry.h"
#include "free_cbspline_chemistry.h"
#include "free_pchip_chemistry.h"
#include "background_chemistry.h"
#include "step_function_chemistry.h"

#include "../config/global_config.h"
#include "../additional/exceptions.h"

#include <vector>
#include <algorithm>
#include <memory>


namespace bear {


//definition of the different chemistry modules with an
//identifier, a keyword to be located in the config file and a short version of the keyword
namespace chemistry_modules{
  enum id {free, iso, eq, cspline, iso_clr, bg, sf, pchip};
  const std::vector<std::string> description {"free", "isoprofile", "equilibrium", "free_cspline", "isoprofile_clr", "background", "step_function", "free_pchip"};
  const std::vector<std::string> description_short {"free", "iso", "eq", "free_cs", "iso_clr", "bg", "sf", "free_p"};
}


inline std::unique_ptr<Chemistry> selectChemistryModule(
  const std::string chemistry_type,
  const std::vector<std::string>& parameters,
  GlobalConfig* config,
  const std::vector<double>& atmos_boundaries)
{
  //find the corresponding chemistry module to the supplied "type" string
  auto it = std::find(
    chemistry_modules::description.begin(),
    chemistry_modules::description.end(),
    chemistry_type);
  auto it_short = std::find(
    chemistry_modules::description_short.begin(),
    chemistry_modules::description_short.end(),
    chemistry_type);


  //no chemistry module is found
  if (it == chemistry_modules::description.end() && it_short == chemistry_modules::description_short.end())
  {
    std::string error_message = "Chemistry type " + chemistry_type + " unknown!\n";
    throw InvalidInput(std::string ("forward_model.config"), error_message);
  }


  //get the id of the chosen module
  chemistry_modules::id module_id = static_cast<chemistry_modules::id>(0);

  if (it != chemistry_modules::description.end())
    module_id = static_cast<chemistry_modules::id>(
      std::distance(chemistry_modules::description.begin(),
      it));
  else
    module_id = static_cast<chemistry_modules::id>(
      std::distance(chemistry_modules::description_short.begin(),
      it_short));


  if (module_id == chemistry_modules::eq)
  {
    if (parameters.size() < 1) {
        std::string error_message = "Equilibrium chemistry requires at least one parameter (the FastChem parameter file)!\n";
        throw InvalidInput(std::string ("forward_model.config"), error_message);}

    const std::vector<std::string> ratio_specs(parameters.begin() + 1, parameters.end());

    return std::make_unique<FastChemChemistry>(
        config->retrieval_folder_path + parameters[0],
        config->nb_omp_processes,
        ratio_specs);
  }

  if (module_id == chemistry_modules::iso)
  {
    return std::make_unique<IsoprofileChemistry>(parameters);
  }

  if (module_id == chemistry_modules::free)
  {
    if (parameters.size() != 3) {
        std::string error_message = "Free chemistry requires exactly three parameters!\n";
        throw InvalidInput(std::string ("forward_model.config"), error_message);}

    return std::make_unique<FreeChemistry>(
        parameters[0],
        std::stoi(parameters[1]),
        std::stoi(parameters[2]),
        atmos_boundaries);
  }

  if (module_id == chemistry_modules::cspline)
  {
    if (parameters.size() != 2) {
        std::string error_message = "Free cubic spline chemistry requires exactly two parameters!\n";
        throw InvalidInput(std::string ("forward_model.config"), error_message);}

    return std::make_unique<FreeCBSplineChemistry>(
        parameters[0],
        std::stoi(parameters[1]));
  }

  if (module_id == chemistry_modules::iso_clr)
  {
    return std::make_unique<IsoprofileCLRChemistry>(parameters);
  }

  if (module_id == chemistry_modules::bg)
  {
    if (parameters.size() != 1) {
        std::string error_message = "Background chemistry requires exactly one parameter!\n";
        throw InvalidInput(std::string ("forward_model.config"), error_message);}

    return std::make_unique<BackgroundChemistry>(parameters[0]);
  }

  if (module_id == chemistry_modules::sf)
  {
    if (parameters.size() != 1) {
        std::string error_message = "Step function chemistry requires exactly one parameter!\n";
        throw InvalidInput(std::string ("forward_model.config"), error_message);}

    return std::make_unique<StepFunctionChemistry>(parameters[0]);
  }

  if (module_id == chemistry_modules::pchip)
  {
    if (parameters.size() != 2) {
        std::string error_message = "Free PCHIP chemistry requires exactly two parameters!\n";
        throw InvalidInput(std::string ("forward_model.config"), error_message);}

    return std::make_unique<FreePchipChemistry>(
        parameters[0],
        std::stoi(parameters[1]));
  }

  //we should never reach this point
  return nullptr;
}


}
#endif
