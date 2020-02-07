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

#include <string>
#include <vector>



namespace fastchem {


template <class double_type>
unsigned int FastChem<double_type>::getChemicalElementIndex(const std::string symbol)
{
  unsigned int index = FASTCHEM_UNKNOWN_SPECIES;

  for (size_t i=0; i<nb_chemical_elements; ++i)
    if (symbol == chemical_elements[i].symbol)
    {
      index = i;
      break;
    }


  return index;
}



template <class double_type>
unsigned int FastChem<double_type>::getMoleculeIndex(const std::string symbol)
{
  unsigned int index = FASTCHEM_UNKNOWN_SPECIES;

  for (size_t i=0; i<molecules.size(); ++i)
    if (symbol == molecules[i].symbol)
    {
      index = i;
      break;
    }


  return index;
}



template <class double_type>
unsigned int FastChem<double_type>::getElementIndex(const std::string symbol)
{
  unsigned int index = FASTCHEM_UNKNOWN_SPECIES;

  for (size_t i=0; i<elements.size(); ++i)
    if (symbol == elements[i].symbol)
    {
      index = i;
      break;
    }


  return index;
}



template <class double_type>
unsigned int FastChem<double_type>::getSpeciesIndex(const std::string symbol)
{
  unsigned int index = FASTCHEM_UNKNOWN_SPECIES;

  for (size_t i=0; i<nb_species; ++i)
    if (symbol == species[i]->symbol)
    {
      index = i;
      break;
    }


  return index;
}



template <class double_type>
std::string FastChem<double_type>::getSpeciesName(const unsigned int species_index)
{
  if (species_index < nb_species)
    return species[species_index]->name;
  else
    return "";
}


template <class double_type>
std::string FastChem<double_type>::getSpeciesSymbol(const unsigned int species_index)
{

  if (species_index < nb_species)
    return species[species_index]->symbol;
  else
    return "";

}



template <class double_type>
std::string FastChem<double_type>::getElementName(const unsigned int species_index)
{
  if (species_index < nb_elements)
    return elements[species_index].name;
  else
    return "";
}



template <class double_type>
std::string FastChem<double_type>::getElementSymbol(const unsigned int species_index)
{

  if (species_index < nb_elements)
    return elements[species_index].symbol;
  else
    return "";

}



template <class double_type>
double FastChem<double_type>::getSpeciesMolecularWeight(const unsigned int species_index)
{

  if (species_index < nb_species)
    return species[species_index]->molecular_weight;
  else
    return 0.;

}



template <class double_type>
double FastChem<double_type>::getElementAbundance(const unsigned int species_index)
{

  if (species_index < nb_elements)
    return elements[species_index].abundance;
  else
    return 0.;

}



template <class double_type>
std::vector<double> FastChem<double_type>::getElementAbundance()
{

  std::vector<double> abundances(nb_elements, 0.0);

  for (size_t i=0; i<nb_elements; ++i)
    abundances[i] = elements[i].abundance;

  return abundances;

}



template class FastChem<double>;
template class FastChem<long double>;

}
