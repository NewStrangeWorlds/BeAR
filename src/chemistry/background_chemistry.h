// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

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
