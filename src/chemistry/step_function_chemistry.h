// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _step_function_chemistry_h
#define _step_function_chemistry_h

#include "chemistry.h"


#include <vector>
#include <string>


namespace bear {



class StepFunctionChemistry : public Chemistry{
  public:
    StepFunctionChemistry(
      const std::string& chemical_species);
    virtual ~StepFunctionChemistry() {}

    virtual bool calcChemicalComposition(
      const std::vector<double>& parameters,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight);
  private:
};


}
#endif 
