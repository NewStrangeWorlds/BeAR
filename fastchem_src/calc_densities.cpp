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
double_type FastChem<double_type>::setInitialHDensity(const double_type total_density, const unsigned int grid_index)
{
  unsigned int H2_ = getMoleculeIndex("H2");
  unsigned int He_ = getElementIndex("He");

  //set initial total total H density
  //we have to treat different cases since H2 or He could be absent
  double_type h_density = 0.0;

  //general case
  if (H2_ != FASTCHEM_UNKNOWN_SPECIES && He_ != FASTCHEM_UNKNOWN_SPECIES)
    h_density = 1./(2. * std::pow(1 + 2.*elements[He_].abundance,2) * std::exp(molecules[H2_].mass_action_constant[grid_index]))
                    * ( (1. + elements[He_].abundance)
                         + 4. * (1 + 2*elements[He_].abundance) * std::exp(molecules[H2_].mass_action_constant[grid_index]) * total_density
                         - std::sqrt(std::pow(1. + elements[He_].abundance, 2)
                         + 4.*(1. + 2.*elements[He_].abundance)*std::exp(molecules[H2_].mass_action_constant[grid_index])*total_density));

  //H2 is present but He is not
  if (H2_ != FASTCHEM_UNKNOWN_SPECIES && He_ == FASTCHEM_UNKNOWN_SPECIES)
    h_density = 1. + 4. * std::exp(molecules[H2_].mass_action_constant[grid_index]) * total_density
                   + std::sqrt( 1. + 4. * std::exp(molecules[H2_].mass_action_constant[grid_index]) * total_density )
                     / (2. * std::exp(molecules[H2_].mass_action_constant[grid_index]));

  //He is present but H2 is not
  if (H2_ == FASTCHEM_UNKNOWN_SPECIES && He_ != FASTCHEM_UNKNOWN_SPECIES)
    h_density = total_density / (1. + elements[He_].abundance);


  //He and H2 are both not present
  if (H2_ == FASTCHEM_UNKNOWN_SPECIES && He_ == FASTCHEM_UNKNOWN_SPECIES)
    h_density = total_density;


  return h_density;
}


template <class double_type>
unsigned int FastChem<double_type>::calcDensity(const double temperature, const double pressure,const unsigned int grid_index,
                                                std::vector<double>& density_n_out, double& h_density_out, double& mean_molecular_weight_out,
                                                std::vector<unsigned int>& element_conserved_out,
                                                unsigned int& nb_pressure_iterations_out, unsigned int& nb_chemistry_iterations_out)
{
  for (auto & i : molecules)  i.calcMassActionConstant(temperature, grid_index);

  //this value will be fixed.
  double_type total_density = pressure/(CONST_K * temperature);

  //initial electron density
  unsigned int e_ = getElementIndex("e-");

  if (e_ != FASTCHEM_UNKNOWN_SPECIES)
    elements[e_].number_density[grid_index] = element_density_minlimit;


  double_type h_density = setInitialHDensity(total_density, grid_index);


  double_type density_iteration_lambda = 0.99; //initial value for lambda
  double_type density_iteration_error = 1.0;   //initial error


  bool fastchem_converged = false;
  bool pressure_converged = false;


  double_type muH = 0;
  double_type amu=1.66055e-24;
  for (size_t i=0; i<nb_elements; ++i)
    muH += elements[i].molecular_weight * chemical_elements[elements[i].element_index].abundance * amu;


  unsigned int nb_iterations = 0;
  unsigned int nb_fastchem_iterations = 0;


  for (nb_iterations=0; nb_iterations<nb_max_pressure_iter; ++nb_iterations)
  {

    fastchem_converged = solveFastchem(temperature, h_density, grid_index, nb_fastchem_iterations);
    pressure_converged = calcTotalHydrogenDensityAlt(temperature, pressure, grid_index,
                                                     h_density, muH, density_iteration_error);
    //pressure_converged = calcTotalHydrogenDensity(temperature, pressure, grid_index,
    //                                              h_density, density_iteration_lambda, density_iteration_error);

    if (std::isnan(h_density)) break;

    if (pressure_converged) break;
  }


  if (!pressure_converged && verbose_level >= 1) std::cout << "Pressure convergence problem in FastChem. :(\n";
  if (!fastchem_converged && verbose_level >= 1) std::cout << "FastChem convergence problem in FastChem. :(\n";


  //return output
  h_density_out = h_density;
  density_n_out.assign(nb_species, 0.0);

  for (size_t i=0; i<nb_species; ++i)
    density_n_out[i] = species[i]->number_density[grid_index];

  mean_molecular_weight_out = calcMeanMolecularWeight(total_density, grid_index);


  for (size_t i=0; i<nb_elements; i++)
    checkElementConservation(elements[i], h_density, grid_index);


  checkChargeConservation(grid_index);



  element_conserved_out.assign(nb_elements, 0);

  for (size_t i=0; i<nb_elements; ++i)
    element_conserved_out[i] = elements[i].element_conserved[grid_index];



  nb_pressure_iterations_out = nb_iterations;
  nb_chemistry_iterations_out = nb_fastchem_iterations;



  unsigned int return_state = FASTCHEM_SUCCESS;

  if (!pressure_converged) return_state = FASTCHEM_NO_PRESSURE_CONVERGENCE;
  if (!fastchem_converged) return_state = FASTCHEM_NO_FASTCHEM_CONVERGENCE;

  if (!fastchem_converged && !pressure_converged) return_state = FASTCHEM_NO_CONVERGENCE;

  return return_state;
}





template <class double_type>
unsigned int FastChem<double_type>::calcDensities(const std::vector<double>& temperature, const std::vector<double>& pressure,
                                                  std::vector < std::vector<double> >& density_out,
                                                  std::vector<double>& h_density_out, std::vector<double>& mean_molecular_weight_out)
{
  std::vector< std::vector<unsigned int> > element_conservation_out;
  std::vector<unsigned int> pressure_iteration_steps_out;
  std::vector<unsigned int> chemistry_iteration_steps_out;
  std::vector<unsigned int> fastchem_flags;


  return calcDensities(temperature, pressure,
                       density_out,
                       h_density_out, mean_molecular_weight_out,
                       element_conservation_out,
                       fastchem_flags,
                       pressure_iteration_steps_out, chemistry_iteration_steps_out);
}




template <class double_type>
unsigned int FastChem<double_type>::calcDensities(const std::vector<double>& temperature, const std::vector<double>& pressure,
                                                  std::vector < std::vector<double> >& density_out,
                                                  std::vector<double>& h_density_out, std::vector<double>& mean_molecular_weight_out,
                                                  std::vector< std::vector<unsigned int> >& element_conserved_out,
                                                  std::vector<unsigned int>& fastchem_flags,
                                                  std::vector<unsigned int>& nb_pressure_iterations_out, std::vector<unsigned int>& nb_chemistry_iterations_out)
{
  if (!is_initialized)
    return FASTCHEM_INITIALIZATION_FAILED;


  unsigned nb_grid_points = temperature.size();


  for (auto & i : species) i->number_density.assign(nb_grid_points, 0);
  for (auto & i : molecules) i.sum.assign(nb_grid_points, 0);
  for (auto & i : molecules) i.mass_action_constant.assign(nb_grid_points, 0);
  for (auto & i : elements) i.element_conserved.assign(nb_grid_points, false);


  element_conserved_out.resize(nb_grid_points);
  nb_pressure_iterations_out.assign(nb_grid_points, 0);
  nb_chemistry_iterations_out.assign(nb_grid_points, 0);


  h_density_out.assign(nb_grid_points, 0);
  mean_molecular_weight_out.assign(nb_grid_points, 0);
  density_out.resize(nb_grid_points);



  std::vector<unsigned int> state(nb_grid_points, 0);


  for (unsigned int i=0; i<nb_grid_points; i++)
    state[i] = calcDensity(temperature[i], pressure[i],
                           i,
                           density_out[i], h_density_out[i], mean_molecular_weight_out[i],
                           element_conserved_out[i],
                           nb_pressure_iterations_out[i], nb_chemistry_iterations_out[i]);


  fastchem_flags = state;


  return *std::max_element(state.begin(),state.end());
}



template <class double_type>
unsigned int FastChem<double_type>::calcDensities(const double temperature, const double pressure,
                                                  std::vector<double>& density_n_out, double& h_density_out, double& mean_molecular_weight_out,
                                                  std::vector<unsigned int>& element_conserved_out,
                                                  unsigned int& nb_pressure_iterations_out, unsigned int& nb_chemistry_iterations_out)
{
  if (!is_initialized)
    return FASTCHEM_INITIALIZATION_FAILED;


  std::vector<double> temperature_vector(1, temperature);
  std::vector<double> pressure_vector(1, pressure);

  std::vector<double> h_density_vector, mean_molecular_weight_vector;
  std::vector< std::vector<double> > density_tensor;

  std::vector< std::vector<unsigned int> > element_conserved_tensor;
  std::vector<unsigned int> nb_pressure_iterations_vector;
  std::vector<unsigned int> nb_chemistry_iterations_vector;
  std::vector<unsigned int> fastchem_flag_vector;



  unsigned int state = calcDensities(temperature_vector, pressure_vector,
                                     density_tensor, h_density_vector, mean_molecular_weight_vector,
                                     element_conserved_tensor,
                                     fastchem_flag_vector, nb_pressure_iterations_vector, nb_chemistry_iterations_vector);

  density_n_out = density_tensor[0];
  h_density_out = h_density_vector[0];
  mean_molecular_weight_out = mean_molecular_weight_vector[0];

  element_conserved_out = element_conserved_tensor[0];
  nb_pressure_iterations_out = nb_pressure_iterations_vector[0];
  nb_chemistry_iterations_out = nb_chemistry_iterations_vector[0];


  return state;
}



template <class double_type>
unsigned int FastChem<double_type>::calcDensities(const double temperature, const double pressure,
                                                  std::vector<double>& density_n_out, double& h_density_out, double& mean_molecular_weight_out)
{
  if (!is_initialized)
    return FASTCHEM_INITIALIZATION_FAILED;


  std::vector<unsigned int> element_conserved;
  unsigned int nb_pressure_iterations, nb_chemistry_iterations;


  unsigned int state = calcDensities(temperature, pressure,
                                     density_n_out, h_density_out, mean_molecular_weight_out,
                                     element_conserved, nb_pressure_iterations, nb_chemistry_iterations);



  return state;
}



template class FastChem<double>;
template class FastChem<long double>;


}


