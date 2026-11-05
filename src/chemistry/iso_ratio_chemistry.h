// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-FileCopyrightText: 2026 Esther Wang
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _iso_ratio_chemistry_h
#define _iso_ratio_chemistry_h

#include "chemistry.h"

#include <vector>
#include <string>


namespace bear {


// Isoprofile chemistry with one isotopologue derived from its parent via a
// retrieved log10 isotopic ratio.
//
// Config syntax:
//   iso_ratio H2O CO 13CO
//
// Parameters (N species → N parameters):
//   params[0..N-2]  linear VMR of each of the first N-1 species
//   params[N-1]     log10 ratio offset from the terrestrial standard
//
// The last listed species is the isotopologue; its abundance is computed as:
//   VMR(iso) = VMR(parent) / standard_ratio * 10^params[N-1]
// where parent is the second-to-last species and standard_ratio = 89.0
// (terrestrial 12C/13C).
class IsoprofileRatioChemistry : public Chemistry {
  public:
    IsoprofileRatioChemistry(const std::vector<std::string>& chemical_species);
    virtual ~IsoprofileRatioChemistry() {}
    virtual bool calcChemicalComposition(
      const std::vector<double>& parameters,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight);

  private:
    static constexpr double standard_ratio = 89.0;  // terrestrial 12C/13C
};


}
#endif
