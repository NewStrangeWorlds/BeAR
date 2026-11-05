// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _fixed_chemistry_h
#define _fixed_chemistry_h

#include "chemistry.h"


#include <vector>
#include <string>


namespace bear {


class FixedChemistry : public Chemistry{
  public:
    FixedChemistry(const std::vector<std::string>& species_symbol);
    virtual ~FixedChemistry() {}
    virtual bool calcChemicalComposition(
      const std::vector<double>& parameters,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight) 
      {
        return true;
      }

    void setChemicalComposition(
      const std::vector<double>& pressure,
      const std::vector<double>& temperature,
      const std::vector<std::vector<double>>& mixing_ratios,
      std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight);
  protected:
    void findSpecies(
      const std::vector<std::string>& species_symbol);
};


}
#endif 