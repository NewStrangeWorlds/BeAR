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

#ifndef _priors_h
#define _priors_h

#include <vector>
#include <string>
#include <iostream>
#include "prior_types.h"


namespace helios {


class Priors{
  public:
    virtual ~Priors();
    void add(
      const std::vector<std::string>& type, 
      const std::vector<std::string>& description, 
      const std::vector<std::vector<double>>& parameter);
    size_t number() {return distributions.size();}
    void printInfo();

    std::vector<BasicPrior*> distributions;
    std::vector<size_t> prior_links;
  protected:
    void addSingle(
      const std::string& type, 
      const std::string& description, 
      const std::vector<double>& parameter);
    void setupLinkedPriors(
      const std::vector<std::string>& type,
      const std::vector<std::string>& description,
      const std::vector<std::vector<double>>& parameter);

  private:
};


}


#endif