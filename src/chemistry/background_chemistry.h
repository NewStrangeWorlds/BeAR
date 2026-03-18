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


#ifndef _background_chemistry_h
#define _background_chemistry_h

#include "chemistry.h"


#include <vector>
#include <string>


namespace bear {


class BackgroundChemistry : public Chemistry{
  public:
    BackgroundChemistry(const std::string& chemical_species);
    virtual ~BackgroundChemistry() {}
    virtual bool calcChemicalComposition(
      const std::vector<double>& parameters,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight);
  protected:
    bool h2he_background = false;
    bool hhe_eq_background = false;

    double mixingRatioBackground(const std::vector<double>& number_densities);

    void backgroundH2He(
      const double mixing_ratio_background, 
      std::vector<double>& number_densities);

    void equilibriumHHe(
      const double pressure, 
      const double temperature,
      const double mixing_ratio_background, 
      std::vector<double>& number_densities);
};


}
#endif 
