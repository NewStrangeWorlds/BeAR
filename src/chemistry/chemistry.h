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


#ifndef _chemistry_h
#define _chemistry_h

#include <vector>


#include "chem_species.h"


namespace bear {


class Chemistry{
  public:
    virtual ~Chemistry() {}
    virtual bool calcChemicalComposition(
      const std::vector<double>& parameters,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight) = 0;
    size_t nbParameters() {return nb_parameters;}
  protected:
    size_t nb_parameters {};
    std::vector<chemical_species_id> species;

    void meanMolecularWeight(
      const std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight) {
        for (size_t i=0; i<number_densities.size(); ++i)
        {
          double mu = 0;

          for (auto & j : constants::species_data)
            mu += number_densities[i][j.id]/number_densities[i][_TOTAL] * j.molecular_weight;

          mean_molecular_weight[i] = mu;
        }
      };

    bool checkMixingRatios(
      const std::vector<std::vector<double>>& number_densities) {
        for (size_t i=0; i<number_densities.size(); ++i)
        {
          double sum = 0;

          for (auto & j : constants::species_data)
            if (j.id != _TOTAL)
              sum += number_densities[i][j.id]/number_densities[i][_TOTAL];

          if (sum > 1+1e-10)
            return false;
        }

        return true;
      };
};


}
#endif 
