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


#ifndef _fastchem_chemistry_h
#define _fastchem_chemistry_h


#include "chemistry.h"

#include "../../_deps/fastchem-src/fastchem_src/fastchem.h"
#include "chem_species.h"

#include <string>
#include <vector>



namespace bear {


class FastChemChemistry : public Chemistry{
  public:
    FastChemChemistry(
      const std::string& fastchen_parameter_file, const size_t nb_openmp_proc);
    virtual ~FastChemChemistry() {}
    
    virtual bool calcChemicalComposition(
      const std::vector<double>& parameters,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight);
  private:
    fastchem::FastChem<long double> fastchem;
    const size_t nb_processes = 0;

    std::vector<double> reference_element_abundances;
    std::vector<size_t> fastchem_species_indices;
};


}
#endif 
