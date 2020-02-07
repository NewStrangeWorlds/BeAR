/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2020 Daniel Kitzmann
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


#ifndef _isoprofile_chemistry_h
#define _isoprofile_chemistry_h

#include "chemistry.h"


#include <vector>
#include <string>


namespace helios {


class IsoprofileChemistry : public Chemistry{
  public:
    IsoprofileChemistry(const std::vector<std::string>& chemical_species);
    virtual ~IsoprofileChemistry() {}
    virtual bool calcChemicalComposition(const std::vector<double>& parameters, const std::vector<double>& temperature, const std::vector<double>& pressure,
                                         std::vector<std::vector<double>>& number_densities, std::vector<double>& mean_molecular_weight);
  protected:
    const double solar_h2 = 0.5;
    const double solar_he = 0.085114;

    bool sodium_free_parameter = false;
};


}
#endif 
