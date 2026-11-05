// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _free_chemistry_h
#define _free_chemistry_h

#include "chemistry.h"

#include "../additional/piecewise_poly.h"


#include <vector>
#include <string>


namespace bear {



class FreeChemistry : public Chemistry{
  public:
    FreeChemistry(
      const std::string& chemical_species,
      const size_t nb_elements_in,
      const size_t polynomial_degree_in,
      const std::vector<double>& atmos_boundaries);
    virtual ~FreeChemistry() {}

    virtual bool calcChemicalComposition(
      const std::vector<double>& parameters,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight);
  private:
    PiecewisePolynomial mixing_ratios;
    const size_t nb_elements {}; 
    const size_t polynomial_degree {};
};


}
#endif 
