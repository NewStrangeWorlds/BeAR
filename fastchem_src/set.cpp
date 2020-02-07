/*
* This file is part of the FastChem code (https://github.com/exoclime/fastchem).
* Copyright (C) 2018 Daniel Kitzmann, Joachim Stock
*
* FastChem is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* FastChem is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* FastChem directory under <license.md>. If not, see
* <http://www.gnu.org/licenses/>.
*/

#include "fastchem.h"

#include <iostream>
#include <string>
#include <vector>



namespace fastchem {


//Set the element abundances for all elements
template <class double_type>
void FastChem<double_type>::setElementAbundance(std::vector<double> abundances)
{

  if (abundances.size() != nb_elements)
  {
    std::cout << "Setting element abundances with an incorrect vector size\n";

    return;
  }


  for (size_t i=0; i<nb_elements; ++i)
  {

    chemical_elements[elements[i].element_index].abundance = abundances[i];
    elements[i].abundance = abundances[i];

  }


  reInitializeFastChem();

}



//Reinitializes certain internal data after the element abundances were changed
template <class double_type>
void FastChem<double_type>::reInitializeFastChem()
{
  element_calculation_order.resize(0);

  //first, we order the elements according to their abundances
  determineElementCalculationOrder();


  //update the definition of the molecular abundances
  for (size_t i=0; i<nb_molecules; ++i)
  {
    molecules[i].abundance = 1e42;

    for (size_t j=0; j<molecules[i].element_indices.size(); ++j)
      if (molecules[i].abundance > elements[molecules[i].element_indices[j]].abundance && elements[molecules[i].element_indices[j]].symbol != "e-")
        molecules[i].abundance = elements[molecules[i].element_indices[j]].abundance;

    //scaled abundances
    molecules[i].abundance_scaled = 1e42;

    for (size_t j=0; j<molecules[i].element_indices.size(); ++j)
    {
      unsigned element_index = molecules[i].element_indices[j];

      if (molecules[i].abundance_scaled > elements[element_index].abundance/molecules[i].stoichometric_vector[element_index]
          && elements[molecules[i].element_indices[j]].symbol != "e-")
        molecules[i].abundance_scaled = elements[element_index].abundance/molecules[i].stoichometric_vector[element_index];
    }
  }


  //update the solver order for the new abundances
  determineSolverOrder();
}




template class FastChem<double>;
template class FastChem<long double>;

}
