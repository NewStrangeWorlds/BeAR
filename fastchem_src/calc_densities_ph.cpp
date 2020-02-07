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


#include <algorithm>
#include <vector>
#include <cmath>


namespace fastchem {


template <class double_type>
unsigned int FastChem<double_type>::calcDensity(const double temperature, const double hydrogen_pressure, const unsigned int grid_index,
                                                std::vector<double>& density_n_out, double& mean_molecular_weight_out,
                                                std::vector<unsigned int>& element_conserved_out,
                                                unsigned int& nb_chemistry_iterations_out)
{
  for (auto & i : molecules)  i.calcMassActionConstant(temperature, grid_index);


  //initial electron density
  unsigned int e_ = getElementIndex("e-");

  if (e_ != FASTCHEM_UNKNOWN_SPECIES)
    elements[e_].number_density[grid_index] = element_density_minlimit;


  double_type h_density = hydrogen_pressure / (CONST_K * temperature);


  bool fastchem_converged = false;


  unsigned int nb_fastchem_iterations = 0;


  fastchem_converged = solveFastchem(temperature, h_density, grid_index, nb_fastchem_iterations);

  if (!fastchem_converged && verbose_level >= 1) std::cout << "FastChem convergence problem in FastChem. :(\n";


  //return output
  density_n_out.assign(nb_species, 0.0);

  for (size_t i=0; i<nb_species; ++i)
    density_n_out[i] = species[i]->number_density[grid_index];


  double_type total_density = 0;

  for (size_t i=0; i<nb_species; ++i)
    total_density += species[i]->number_density[grid_index];

  mean_molecular_weight_out = calcMeanMolecularWeight(total_density, grid_index);


  for (size_t i=0; i<nb_elements; i++)
    checkElementConservation(elements[i], h_density, grid_index);


  checkChargeConservation(grid_index);



  element_conserved_out.assign(nb_elements, 0);

  for (size_t i=0; i<nb_elements; ++i)
    element_conserved_out[i] = elements[i].element_conserved[grid_index];

  nb_chemistry_iterations_out = nb_fastchem_iterations;



  unsigned int return_state = FASTCHEM_SUCCESS;

  if (!fastchem_converged) return_state = FASTCHEM_NO_FASTCHEM_CONVERGENCE;


  return return_state;
}




template <class double_type>
unsigned int FastChem<double_type>::calcDensities(const std::vector<double>& temperature, const std::vector<double>& hydrogen_pressure,
                                                  std::vector < std::vector<double> >& density_out,
                                                  std::vector<double>& mean_molecular_weight_out,
                                                  std::vector< std::vector<unsigned int> >& element_conserved_out,
                                                  std::vector<unsigned int>& fastchem_flags,
                                                  std::vector<unsigned int>& nb_chemistry_iterations_out)
{
  if (!is_initialized)
    return FASTCHEM_INITIALIZATION_FAILED;


  unsigned nb_grid_points = temperature.size();


  for (auto & i : species) i->number_density.assign(nb_grid_points, 0);
  for (auto & i : molecules) i.sum.assign(nb_grid_points, 0);
  for (auto & i : molecules) i.mass_action_constant.assign(nb_grid_points, 0);
  for (auto & i : elements) i.element_conserved.assign(nb_grid_points, false);


  element_conserved_out.resize(nb_grid_points);
  nb_chemistry_iterations_out.assign(nb_grid_points, 0);


  mean_molecular_weight_out.assign(nb_grid_points, 0);
  density_out.resize(nb_grid_points);



  std::vector<unsigned int> state(nb_grid_points, 0);


  for (unsigned int i=0; i<nb_grid_points; i++)
    state[i] = calcDensity(temperature[i], hydrogen_pressure[i],
                           i,
                           density_out[i], mean_molecular_weight_out[i],
                           element_conserved_out[i],
                           nb_chemistry_iterations_out[i]);


  fastchem_flags = state;


  return *std::max_element(state.begin(),state.end());
}



template class FastChem<double>;
template class FastChem<long double>;


}


