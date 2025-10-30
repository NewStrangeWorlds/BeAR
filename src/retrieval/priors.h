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

#ifndef _priors_h
#define _priors_h

#include <vector>
#include <string>
#include <iostream>
#include "prior_types.h"


namespace bear {


struct PriorConfig{
  PriorConfig(
    const std::string& type_, 
    const std::string& description_, 
    const std::vector<double>& parameter_,
    const std::string& unit_)
    : type(type_), description(description_), parameter(parameter_), unit(unit_) {}
  PriorConfig(
    const std::string& type_, 
    const std::string& description_, 
    const std::vector<double>& parameter_)
    : type(type_), description(description_), parameter(parameter_) {}
  
  std::string type = "";
  std::string description = "";
  std::vector<double> parameter;
  std::string unit = "";
};



class Priors{
  public:
    virtual ~Priors();

    void init(
      const std::string& folder_path, 
      const size_t nb_total_param);
    void init(
      const std::vector<PriorConfig>& priors_config,
      const size_t nb_total_param);

    void add(
      const std::vector<PriorConfig>& priors_config);
    
    size_t number() {return distributions.size();}
    void printInfo();

    std::vector<BasicPrior*> distributions;
    std::vector<size_t> prior_links;
  protected:
    void readConfigFile(
      const std::string& file_path, 
      std::vector<std::string>& prior_type, 
      std::vector<std::string>& prior_description, 
      std::vector<std::vector<double>>& prior_parameter,
      std::vector<std::string>& prior_unit);
    void addSingle(
      const std::string& type, 
      const std::string& description, 
      const std::vector<double>& parameter,
      const std::string& unit);
    void setupLinkedPriors(
      const std::vector<std::string>& type,
      const std::vector<std::string>& description,
      const std::vector<std::vector<double>>& parameter);
};


}


#endif