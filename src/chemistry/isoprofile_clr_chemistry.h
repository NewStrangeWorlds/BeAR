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


#ifndef _isoprofile_clr_chemistry_h
#define _isoprofile_clr_chemistry_h

#include "chemistry.h"


#include <vector>
#include <string>


namespace helios {


//isoprofiles with a centred-log-ratio prior
//does not use a buffer gas like the normal isoprofile chemistry
class IsoprofileCLRChemistry : public Chemistry{
  public:
    IsoprofileCLRChemistry(const std::vector<std::string>& chemical_species);
    virtual ~IsoprofileCLRChemistry() {}
    virtual bool calcChemicalComposition(
      const std::vector<double>& parameters,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight);
  protected:
    const double min_number_density = 1e-12;
    bool transformPriors(
      const std::vector<double>& priors,
      std::vector<double>& mixing_ratios);
};


}
#endif 
