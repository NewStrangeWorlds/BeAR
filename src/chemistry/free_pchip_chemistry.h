// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _free_pchip_chemistry_h
#define _free_pchip_chemistry_h

#include "chemistry.h"

#include <vector>
#include <string>


namespace bear {


class FreePchipChemistry : public Chemistry{
  public:
    FreePchipChemistry(
      const std::string& chemical_species,
      const size_t nb_control_points_);
    virtual ~FreePchipChemistry() {}

    virtual bool calcChemicalComposition(
      const std::vector<double>& parameters,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight);
  private:
    const size_t nb_control_points = 0;
};


}
#endif
