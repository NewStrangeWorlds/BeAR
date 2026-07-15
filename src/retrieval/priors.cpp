/*
* This file is part of the BeAR code (https://github.com/newstrangeworlds/BeAR).
* Copyright (C) 2025 Daniel Kitzmann
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
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <sstream>

#include "priors.h"

#include "prior_types.h"
#include "../additional/exceptions.h"


namespace bear{


Priors::~Priors()
{
   for (size_t i=0; i<distributions.size(); ++i)
    delete distributions[i];
}



void Priors::init(
  const std::string& folder_path, 
  const size_t nb_total_param)
{
  const std::string file_name = folder_path + "priors.config";

  std::vector<std::string> prior_type; 
  std::vector<std::string> prior_description; 
  std::vector<std::vector<double>> prior_parameter;
  std::vector<std::string> prior_unit;

  readConfigFile(file_name, prior_type, prior_description, prior_parameter, prior_unit);

  if (prior_type.size() != nb_total_param)
  {
    std::string error_message = "Found " 
      + std::to_string(prior_type.size()) 
      + " priors in priors.config but expected " 
      + std::to_string(nb_total_param) + "\n";
    
    throw InvalidInput(std::string ("Priors::init"), error_message);
  }
  
  std::vector<PriorConfig> priors_config;
  
  for (size_t i=0; i<prior_type.size(); ++i)
    priors_config.push_back(PriorConfig(prior_type[i], prior_description[i], prior_parameter[i], prior_unit[i]));

  add(priors_config);
}



void Priors::init(
  const std::vector<PriorConfig>& priors_config,
  const size_t nb_total_param)
{
  if (priors_config.size() != nb_total_param)
  {
    std::string error_message = "Found " 
      + std::to_string(priors_config.size()) 
      + " priors config but expected " 
      + std::to_string(nb_total_param) + "\n";
    
    throw InvalidInput(std::string ("Priors::init"), error_message);
  }

  add(priors_config);
}



void Priors::printInfo()
{
  std::cout << "\n" << "List of priors: \n";

  for (auto & i : distributions)
    i->printInfo();

  std::cout << "\n";
}



//Writes a small legend that maps each column of the posterior file to the prior
//it corresponds to. The posterior rows are laid out as
//  p_0  p_1  ...  p_(numberFree()-1)  log_likelihood
//with the parameter columns in the order of the compressed sampling cube.
void Priors::writeParameterList(const std::string& file_path)
{
  std::fstream file(file_path.c_str(), std::ios::out);

  if (file.fail())
  {
    std::cout << "Warning: could not write parameter list to " << file_path << "\n";
    return;
  }

  const size_t nb_free = numberFree();

  file << "# BeAR posterior parameter legend\n";
  file << "# Each sampled parameter below is one column of the posterior file,\n";
  file << "# in this order:  p_0  p_1  ...  p_" << (nb_free == 0 ? 0 : nb_free - 1)
       << "  log_likelihood\n";
  file << "# (the final posterior column is the log-likelihood and is not listed here)\n";
  file << "#\n";
  file << "# column\tparameter\n";

  size_t column = 0;

  for (size_t i=0; i<distributions.size(); ++i)
  {
    //fixed (delta) priors have free_cube_index < 0; linked priors share their
    //target's column and do not add one of their own
    if (free_cube_index[i] >= 0
        && distributions[i]->distributionType() != "Linked prior")
    {
      file << column << "\t" << distributions[i]->parameterName() << "\n";
      ++column;
    }
  }

  //for completeness, note the priors that do not get their own posterior column
  bool has_non_sampled = false;

  for (size_t i=0; i<distributions.size(); ++i)
    if (distributions[i]->isFixed()
        || distributions[i]->distributionType() == "Linked prior")
      has_non_sampled = true;

  if (has_non_sampled)
  {
    file << "#\n# not sampled (no posterior column of their own):\n";

    for (size_t i=0; i<distributions.size(); ++i)
    {
      if (distributions[i]->distributionType() == "Linked prior")
        file << "#   " << distributions[i]->parameterName()
             << "\tlinked (shares column " << free_cube_index[i] << ")\n";
      else if (distributions[i]->isFixed())
        file << "#   " << distributions[i]->parameterName() << "\tfixed\n";
    }
  }

  file.close();

  std::cout << "Wrote posterior parameter legend to " << file_path << "\n";
}



void Priors::add(
  const std::vector<PriorConfig>& priors_config)
{
  std::vector<std::string> type(priors_config.size(), "");
  std::vector<std::string> description(priors_config.size(), "");
  std::vector<std::vector<double>> parameter(priors_config.size(), std::vector<double>());
  std::vector<std::string> unit(priors_config.size(), "");
  std::vector<std::string> link_target(priors_config.size(), "");

  for (size_t i=0; i<priors_config.size(); ++i)
  {
    type[i] = priors_config[i].type;
    description[i] = priors_config[i].description;

    for (size_t j=0; j<priors_config[i].parameter.size(); ++j)
      parameter[i].push_back(priors_config[i].parameter[j]);

    unit[i] = priors_config[i].unit;
    link_target[i] = priors_config[i].link_target;
  }

  for (size_t i=0; i<priors_config.size(); ++i)
    addSingle(type[i], description[i], parameter[i], unit[i]);

  setupLinkedPriors(type, description, link_target);
}



void Priors::addSingle(
  const std::string& type, 
  const std::string& description, 
  const std::vector<double>& parameter,
  const std::string& unit)
{
  if (type == "uniform")
  {
    if (parameter.size() != 2)
    {
      std::string error_message = "uniform prior " + description + " requires two parameters with an optional unit!\n";
      throw InvalidInput(std::string ("pirors.config"), error_message);
    }

    if (unit == "")
    {
      UniformPrior* uniform_prior = new UniformPrior(
        description, parameter[0], parameter[1], "");
      distributions.push_back(uniform_prior);
    }
    else
    {
      UniformPrior* uniform_prior = new UniformPrior(
        description, parameter[0], parameter[1], unit);
      distributions.push_back(uniform_prior);
    }

    return;
  }


  if (type == "log_uniform")
  {
    if (parameter.size() != 2)
    {
      std::string error_message = "log uniform prior " + description + " requires two parameters with an optional unit!\n";
      throw InvalidInput(std::string ("pirors.config"), error_message);
    }

    if (unit == "")
    {
      LogUniformPrior* uniform_prior = new LogUniformPrior(
        description, parameter[0], parameter[1], "");
      distributions.push_back(uniform_prior);
    }
    else
    {
      LogUniformPrior* uniform_prior = new LogUniformPrior(
        description, parameter[0], parameter[1], unit);
      distributions.push_back(uniform_prior);
    }

    return;
  }


  if (type == "gaussian")
  {
    if (parameter.size() != 2)
    {
      std::string error_message = "Gaussian prior " + description + " requires two parameters with an optional unit!\n";
      throw InvalidInput(std::string ("pirors.config"), error_message);
    }

    if (unit == "")
    {
      GaussianPrior* gaussian_prior = new GaussianPrior(
        description, parameter[0], parameter[1], "");
      distributions.push_back(gaussian_prior);
    }
    else
    {
      GaussianPrior* gaussian_prior = new GaussianPrior(
        description, parameter[0], parameter[1], unit);
      distributions.push_back(gaussian_prior);
    }

    return;
  }


  if (type == "delta")
  {
    if (parameter.size() != 1)
    {
      std::string error_message = "Delta distribution prior " + description + " requires one parameter with an optional unit!\n";
      throw InvalidInput(std::string ("pirors.config"), error_message);
    }

    if (unit == "")
    {
      DeltaPrior* delta_prior = new DeltaPrior(
        description, parameter[0], "");
      distributions.push_back(delta_prior);
    }
    else
    {
      DeltaPrior* delta_prior = new DeltaPrior(
        description, parameter[0], unit);
      distributions.push_back(delta_prior);
    }

    return;
  }


  if (type == "linked")
  {
    //Placeholder only; the real LinkedPrior is created in setupLinkedPriors once
    //every distribution exists and the target name can be resolved. The value
    //here is irrelevant since the placeholder is deleted there.
    DeltaPrior* delta_prior = new DeltaPrior(description, 0.0, "");
    distributions.push_back(delta_prior);

    return;
  }

  std::string error_message = "Prior type " + type + " for " + description + " unknown!\n";
  throw InvalidInput(std::string ("Retrieval::setPrior"), error_message);
}




void Priors::setupLinkedPriors(
  const std::vector<std::string>& type,
  const std::vector<std::string>& description,
  const std::vector<std::string>& link_target)
{
  //add() may run more than once (model parameters, then the error-inflation
  //prior). Grow prior_links to cover all distributions while preserving earlier
  //entries — assigning to type.size() here would truncate a previous batch's
  //links and corrupt the free_cube_index lookup below. Linked priors only ever
  //occur in the first batch, where the batch index equals the global index.
  prior_links.resize(distributions.size(), 0);

  for (size_t i=0; i<type.size(); ++i)
  {
    if (type[i] == "linked")
    {
      //resolve the link target by NAME to its position in the canonical order
      size_t prior_index = 0;
      bool found = false;

      for (size_t j=0; j<description.size(); ++j)
        if (description[j] == link_target[i])
        {
          prior_index = j;
          found = true;
          break;
        }

      if (!found)
      {
        std::string error_message = "Linked prior '" + description[i]
          + "' references unknown target parameter '" + link_target[i] + "'!\n";
        throw InvalidInput(std::string ("priors.config"), error_message);
      }

      if (prior_index == i)
      {
        std::string error_message = "Linked prior '" + description[i] + "' cannot link to itself!\n";
        throw InvalidInput(std::string ("priors.config"), error_message);
      }

      if (type[prior_index] == "linked")
      {
        std::string error_message = "Linked prior '" + description[i] + "' cannot link to another linked prior!\n";
        throw InvalidInput(std::string ("priors.config"), error_message);
      }

      //delete the place holder now that the target is validated
      delete(distributions[i]);

      LinkedPrior* linked_prior = new LinkedPrior(description[i], distributions[prior_index]);
      distributions[i] = linked_prior;
      prior_links[i] = prior_index;
    }
  }

  // Compute free_cube_index: maps each prior to its position in the compressed
  // free-parameter cube, or -1 if the prior is fixed (delta).
  free_cube_index.assign(distributions.size(), -1);
  int j = 0;
  for (size_t i = 0; i < distributions.size(); i++) {
    if (!distributions[i]->isFixed() &&
        distributions[i]->distributionType() != "Linked prior")
      free_cube_index[i] = j++;
  }
  for (size_t i = 0; i < distributions.size(); i++) {
    if (distributions[i]->distributionType() == "Linked prior")
      free_cube_index[i] = free_cube_index[prior_links[i]];
  }
}


size_t Priors::numberFree() const {
  //Number of distinct slots in the compressed free-parameter cube. Linked priors
  //share their target's slot, so counting every free_cube_index >= 0 would
  //over-count them as extra dimensions; take (max index + 1) instead. With no
  //linked priors this is identical to counting the non-fixed distributions.
  int max_index = -1;
  for (int idx : free_cube_index)
    if (idx > max_index) max_index = idx;
  return static_cast<size_t>(max_index + 1);
}


std::vector<double> Priors::expandFreeToFull(const std::vector<double>& free_phys) const {
  std::vector<double> full(distributions.size(), 0.0);
  for (size_t i = 0; i < distributions.size(); i++) {
    if (free_cube_index[i] < 0)
      full[i] = distributions[i]->parameterPhysicalValue(0.0);
    else
      full[i] = free_phys[free_cube_index[i]];
  }
  return full;
}



//Parse priors.config into a name-keyed map. Unlike readConfigFile this is
//order-independent: the "description" column is used as the parameter name/key.
//Blank lines and lines starting with '#' are skipped.
std::map<std::string, PriorConfig> Priors::parseConfigToMap(
  const std::string& file_path)
{
  std::fstream file;
  file.open(file_path.c_str(), std::ios::in);

  if (file.fail())
    throw FileNotFound(std::string ("Priors::parseConfigToMap"), file_path);

  auto is_number = [](const std::string& s){
    std::istringstream iss(s);
    double d;
    return iss >> std::noskipws >> d && iss.eof();};

  std::map<std::string, PriorConfig> prior_map;
  std::string line;

  while (std::getline(file, line))
  {
    std::istringstream input(line);

    std::string type, description;
    input >> type >> description;

    if (type.empty() || type[0] == '#')
      continue;

    if (prior_map.count(description))
    {
      std::string error_message = "Duplicate prior name '" + description + "' in priors.config!\n";
      throw InvalidInput(std::string ("priors.config"), error_message);
    }

    //A linked prior carries the NAME of the parameter it links to (not numbers).
    //Resolution to a slot happens later against the canonical parameter order.
    if (type == "linked")
    {
      std::string target;
      input >> target;

      if (target.empty())
      {
        std::string error_message = "linked prior '" + description
          + "' requires a target parameter name!\n";
        throw InvalidInput(std::string ("priors.config"), error_message);
      }

      PriorConfig config(type, description, std::vector<double>{}, "");
      config.link_target = target;
      prior_map.emplace(description, config);

      continue;
    }

    std::vector<std::string> parameter;
    std::string single_parameter;

    while (input >> single_parameter)
      parameter.push_back(single_parameter);

    if (parameter.empty())
    {
      std::string error_message = "Prior '" + description + "' has no parameters!\n";
      throw InvalidInput(std::string ("priors.config"), error_message);
    }

    std::vector<double> double_parameter;
    std::string unit = "";

    const size_t nb_numbers = is_number(parameter.back())
      ? parameter.size() : parameter.size() - 1;

    for (size_t i=0; i<nb_numbers; ++i)
    {
      if (is_number(parameter[i]) == false)
      {
        std::string error_message = "Prior value " + parameter[i] + " cannot be converted into a number!\n";
        throw InvalidInput(std::string ("priors.config"), error_message);
      }

      double_parameter.push_back(std::stod(parameter[i]));
    }

    if (nb_numbers < parameter.size())
      unit = parameter.back();

    prior_map.emplace(
      description,
      PriorConfig(type, description, double_parameter, unit));
  }

  file.close();

  return prior_map;
}



//Build the prior distributions in the canonical order given by ordered_names,
//looking each up by name in the parsed map. Missing or stray priors are hard
//errors with a descriptive message - replacing the old positional count check.
void Priors::initFromMap(
  const std::map<std::string, PriorConfig>& prior_map,
  const std::vector<std::string>& ordered_names)
{
  std::vector<PriorConfig> priors_config;
  priors_config.reserve(ordered_names.size());

  for (const auto& name : ordered_names)
  {
    auto it = prior_map.find(name);

    if (it == prior_map.end())
    {
      std::string error_message =
        "No prior found in priors.config for required parameter '" + name + "'.\n"
        "Expected parameters (in canonical order):\n";

      for (const auto& n : ordered_names)
        error_message += "  " + n + "\n";

      throw InvalidInput(std::string ("Priors::initFromMap"), error_message);
    }

    priors_config.push_back(it->second);
  }

  //flag any priors in the file that don't correspond to a model parameter
  //(typo guard - this was silently impossible to detect with positional parsing)
  const std::set<std::string> expected(ordered_names.begin(), ordered_names.end());

  for (const auto& entry : prior_map)
    if (expected.count(entry.first) == 0)
    {
      std::string error_message =
        "Prior '" + entry.first + "' in priors.config matches no model parameter.\n";

      throw InvalidInput(std::string ("Priors::initFromMap"), error_message);
    }

  add(priors_config);
}



void Priors::readConfigFile(
  const std::string& file_path, 
  std::vector<std::string>& prior_type, 
  std::vector<std::string>& prior_description, 
  std::vector<std::vector<double>>& prior_parameter,
  std::vector<std::string>& prior_unit)
{

  std::fstream file;
  file.open(file_path.c_str(), std::ios::in);

  if (file.fail())  
    throw FileNotFound(std::string ("ForwardModel::readPriorConfigFile"), file_path);


  auto is_number = [](const std::string& s){
    std::istringstream iss(s);
    double d;
    return iss >> std::noskipws >> d && iss.eof();};


  std::string line;

  while (std::getline(file, line))
  {
    std::istringstream input(line);

    std::string type, description;
    std::vector<std::string> parameter;
    std::string unit = "";

    input >> type >> description;

    std::string single_parameter;

    while (input >> single_parameter)
      parameter.push_back(single_parameter);

    std::vector<double> double_parameter;

    if (is_number(parameter.back()) == true)
    {
      for (size_t i=0; i<parameter.size(); ++i)
      {
        if (is_number(parameter[i]) == false)
        {
          std::string error_message = "Prior value " + parameter[i] + " cannot be converted into a number!\n";
          throw InvalidInput(std::string ("pirors.config"), error_message);
        }
        
        double_parameter.push_back(std::stod(parameter[i]));
      }
        
    }
    else
    {
      for (size_t i=0; i<parameter.size()-1; ++i)
      {
        if (is_number(parameter[i]) == false)
        {
          std::string error_message = "Prior value " + parameter[i] + " cannot be converted into a number!\n";
          throw InvalidInput(std::string ("pirors.config"), error_message);
        }

        double_parameter.push_back(std::stod(parameter[i]));
      }
      
      unit = parameter.back();
    }

    prior_type.push_back(type);
    prior_description.push_back(description);
    prior_parameter.push_back(double_parameter);
    prior_unit.push_back(unit);
  }

  file.close();
}


}
