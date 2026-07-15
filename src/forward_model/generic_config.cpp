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
#include <sstream>
#include <iomanip>
#include <limits>
#include <vector>
#include <algorithm>

#include "generic_config.h"
#include "../additional/exceptions.h"


namespace bear{


toml::table GenericConfig::parseConfigFile(const std::string& file_name)
{
  try
  {
    return toml::parse_file(file_name);
  }
  catch (const toml::parse_error& e)
  {
    std::ostringstream error_message;
    error_message << "Error parsing " << file_name << ": "
                  << e.description() << " (at " << e.source().begin << ")\n";
    throw InvalidInput(std::string ("GenericConfig::readConfigFile"), error_message.str());
  }
}



std::string GenericConfig::nodeToString(const toml::node& node)
{
  if (auto v = node.value<std::string>())
    return *v;

  if (auto v = node.value<int64_t>())
    return std::to_string(*v);

  if (auto v = node.value<double>())
  {
    std::ostringstream os;
    os << std::setprecision(std::numeric_limits<double>::max_digits10) << *v;
    return os.str();
  }

  if (auto v = node.value<bool>())
    return *v ? "true" : "false";

  return "";
}



std::string GenericConfig::missingSectionMessage(const std::string& key)
{
  return "Required section '" + key + "' is missing or empty in the forward model config.\n"
    "If it is written after a [section] header, that header captures it: in TOML a\n"
    "[section] absorbs every following key until the next header. Put bare inline\n"
    "keys (e.g. " + key + " = [...]) before any [section] header, or give this\n"
    "section its own header too.\n";
}



double GenericConfig::requirePositive(
  const toml::table& tbl,
  const std::string& section,
  const std::string& key)
{
  auto value = tbl[key].value<double>();

  if (!value)
    throw InvalidInput(std::string ("GenericConfig::requirePositive"),
      "Key '" + key + "' in [" + section + "] is missing or not a number "
      "(check for a typo in the key name).\n");

  if (*value <= 0.0)
    throw InvalidInput(std::string ("GenericConfig::requirePositive"),
      "Key '" + key + "' in [" + section + "] must be a positive number.\n");

  return *value;
}



void GenericConfig::readAtmosphereConfig(
  const toml::table& cfg,
  size_t& nb_grid_points,
  std::vector<double>& pressure_boundaries)
{
  const toml::table* tbl = cfg["atmosphere"].as_table();

  if (!tbl)
    throw InvalidInput(std::string ("GenericConfig::readAtmosphereConfig"),
      std::string("Missing '[atmosphere]' table in config file\n"));

  nb_grid_points = static_cast<size_t>(requirePositive(*tbl, "atmosphere", "nb_grid_points"));

  pressure_boundaries = {
    requirePositive(*tbl, "atmosphere", "bottom_pressure"),
    requirePositive(*tbl, "atmosphere", "top_pressure")};

  std::cout << "- Atmosphere levels: " << nb_grid_points << "\n";
  std::cout << "- Pressure boundaries: "
            << pressure_boundaries[0] << " " << pressure_boundaries[1] << "\n";
}



void GenericConfig::readModelBlock(
  const toml::table& cfg,
  const std::string& key,
  std::string& model,
  std::vector<std::string>& parameters,
  const std::string& log_label)
{
  const toml::table* tbl = cfg[key].as_table();

  if (!tbl)
    throw InvalidInput(std::string ("GenericConfig::readModelBlock"),
      "Missing or invalid '[" + key + "]' table in config file\n");

  model = (*tbl)["model"].value_or(std::string(""));

  if (const toml::array* params = (*tbl)["params"].as_array())
    for (auto&& p : *params)
      parameters.push_back(nodeToString(p));

  std::cout << "- " << log_label << ": " << model;
  for (auto & i : parameters) std::cout << "  " << i;
  std::cout << "\n";
}



void GenericConfig::readModelList(
  const toml::table& cfg,
  const std::string& key,
  std::vector<std::string>& models,
  std::vector<std::vector<std::string>>& parameters,
  bool skip_none,
  bool required,
  const std::string& log_label)
{
  //absent or empty array means "none" (allowed unless the section is required)
  if (const toml::array* arr = cfg[key].as_array())
  {
    for (auto&& entry : *arr)
    {
      const toml::table* t = entry.as_table();
      if (!t) continue;

      std::string model = (*t)["model"].value_or(std::string(""));

      std::cout << "- " << log_label << ": " << model << "\n";

      if (skip_none &&
          (model == "None" || model == "none" || model == "No"
           || model == "no" || model == "N" || model == "n"))
        continue;

      models.push_back(model);
      parameters.resize(parameters.size() + 1);

      if (const toml::array* params = (*t)["params"].as_array())
        for (auto&& p : *params)
          parameters.back().push_back(nodeToString(p));
    }
  }

  if (required && models.empty())
    throw InvalidInput(std::string ("GenericConfig::readModelList"),
      missingSectionMessage(key));
}



void GenericConfig::readOpacityConfig(
  const toml::table& cfg,
  const std::string& key,
  std::vector<std::string>& species_symbol,
  std::vector<std::string>& species_folder,
  bool required)
{
  if (const toml::array* arr = cfg[key].as_array())
  {
    for (auto&& entry : *arr)
    {
      const toml::table* t = entry.as_table();
      if (!t) continue;

      std::string species = (*t)["species"].value_or(std::string(""));
      std::string folder  = (*t)["folder"].value_or(std::string(""));

      if (!species.empty() && !folder.empty())
      {
        species_symbol.push_back(species);
        species_folder.push_back(folder);
      }
    }
  }

  if (required && species_symbol.empty())
    throw InvalidInput(std::string ("GenericConfig::readOpacityConfig"),
      missingSectionMessage(key));

  std::cout << "- Opacity species:\n";
  for (size_t i=0; i<species_symbol.size(); ++i)
    std::cout << "   species " << species_symbol[i] << "\t folder: " << species_folder[i] << "\n";

  std::cout << "\n";
}



bool GenericConfig::readBooleanParameter(
  const toml::table& cfg,
  const std::string& key,
  bool default_value)
{
  auto node = cfg[key];

  bool value = default_value;

  if (auto b = node.value<bool>())
    value = *b;
  else if (auto s = node.value<std::string>())
  {
    //tolerate legacy Yes/No spellings in hand-written TOML
    if (*s == "Yes" || *s == "yes" || *s == "Y" || *s == "y")
      value = true;
    else if (*s == "No" || *s == "no" || *s == "N" || *s == "n")
      value = false;
    else
      throw InvalidInput(std::string ("GenericConfig::readBooleanParameter"),
        "Boolean parameter value for " + key + " in config file is invalid\n");
  }
  else if (node)
    throw InvalidInput(std::string ("GenericConfig::readBooleanParameter"),
      "Boolean parameter value for " + key + " in config file is invalid\n");
  //absent -> keep default_value

  std::cout << "  - " << key << ": " << (value ? "yes" : "no") << "\n";

  return value;
}



std::string GenericConfig::readParameter(
  const toml::table& cfg,
  const std::string& key,
  const std::vector<std::string>& allowed_values)
{
  std::string param = cfg[key].value_or(std::string(""));

  auto it = std::find(allowed_values.begin(), allowed_values.end(), param);

  if (it == allowed_values.end())
    throw InvalidInput(std::string ("GenericConfig::readParameter"),
      "Parameter value " + param + " for " + key + " in config file is invalid\n");

  return param;
}



std::vector<chemical_species_id> GenericConfig::readChemicalSpecies(
  const toml::table& cfg,
  const std::string& key)
{
  std::vector<chemical_species_id> species_to_save;

  if (const toml::array* arr = cfg[key].as_array())
  {
    for (auto&& entry : *arr)
    {
      const std::string species = nodeToString(entry);

      for (size_t j=0; j<constants::species_data.size(); ++j)
      {
        if (constants::species_data[j].symbol == species)
        {
          species_to_save.push_back(constants::species_data[j].id);
          break;
        }
      }
    }
  }

  std::cout << "  - " << key << ": ";
  for (auto & i : species_to_save)
    std::cout << constants::species_data[i].symbol << " ";
  std::cout << "\n";

  return species_to_save;
}


}
