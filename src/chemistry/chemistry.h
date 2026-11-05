// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _chemistry_h
#define _chemistry_h

#include <vector>
#include <string>


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
    size_t nbParameters() {return parameter_names.size();}
    const std::vector<std::string>& parameterNames() const {return parameter_names;}
  protected:
    std::vector<std::string> parameter_names;
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

          if (sum > 1.1+1e-10)
            return false;
        }

        return true;
      };
};


}
#endif 
