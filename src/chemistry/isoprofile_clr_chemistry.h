// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _isoprofile_clr_chemistry_h
#define _isoprofile_clr_chemistry_h

#include "chemistry.h"


#include <vector>
#include <string>


namespace bear {


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
