/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2019 Daniel Kitzmann
*
* Helios-r2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Helios-r2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* Helios-r2 directory under <LICENSE>. If not, see
* <http://www.gnu.org/licenses/>.
*/


#include "../forward_model/piecewise_poly.h"


#include <cstddef>
#include <vector>
#include <iostream>
#include <cmath>


namespace helios{

//Gauss-Lobatto quadrature nodes on the reference [-1, +1] element
//for orders up to six
const std::vector<double> gl_0 = {0.0};
const std::vector<double> gl_1 = {-1.0, 1.0};
const std::vector<double> gl_2 = {-1.0, 0.0, 1.0};
const std::vector<double> gl_3 = {-1.0, -0.447214, 0.447214, 1.0};
const std::vector<double> gl_4 = {-1.0, -0.654654, 0.654654, 1.0};
const std::vector<double> gl_5 = {-1.0, -0.654654, 0.0, 0.654654, 1.0};
const std::vector<double> gl_6 = {-1.0, -7.65055 -0.285232, 0.285232, 0.765055, 1.0};

const std::vector< std::vector<double> > quadrature_nodes{gl_0, gl_1, gl_2, gl_3, gl_4, gl_5, gl_6};



  
PiecewisePolynomial::PiecewisePolynomial(const size_t element_number, const size_t polynomial_order, const double domain_boundaries [2])
{
  double log_boundaries [2];
  log_boundaries[0] = std::log10(domain_boundaries[0]);
  log_boundaries[1] = std::log10(domain_boundaries[1]);

  nb_elements = element_number;
  if (nb_elements < 1) nb_elements = 1;

  nb_edges = nb_elements + 1;


  order = polynomial_order;
  if (polynomial_order < 1) order = 1;
  if (polynomial_order > 6) order = 6;


  createElementGrid(log_boundaries);  
}



size_t PiecewisePolynomial::createElementGrid(const double domain_boundaries [2])
{

  double domain_size = domain_boundaries[0] - domain_boundaries[1];
  double element_size = domain_size / nb_elements;

 
  element_edges.resize(nb_elements+1);

  element_edges.front() = domain_boundaries[0];
  element_edges.back() = domain_boundaries[1];

  
  for (size_t i=1; i<nb_edges-1; ++i)
    element_edges[i] = element_edges[i-1] - element_size;

  
  elements.reserve(nb_elements);

  for (size_t i=0; i<nb_elements; ++i)
  {
    double edges [2] = {element_edges[i], element_edges[i+1]};

    Element element(edges, order);

    elements.push_back(element);
  }


  for (size_t i=0; i<nb_elements; ++i)
    for (size_t j=0; j<elements[i].nb_dof-1; ++j)
      dof_vertices.push_back(elements[i].dof_vertices[j]);

  dof_vertices.push_back(elements.back().dof_vertices.back());

  nb_dof = dof_vertices.size();
  dof_values.assign(nb_dof, 0.0);


  return 1;
}




void PiecewisePolynomial::setDOFvalues(const std::vector<double>& values)
{
  
  if (values.size() != nb_dof)
  {
     std::cout << "Passed vector length does not correspond to the number of dof!\n";
     return;
  }

  dof_values = values;

  //set the dof values in each element
  size_t global_dof_index = 0;

  for (size_t i=0; i<nb_elements; ++i)
  {
    
    for (size_t j=0; j<elements[i].nb_dof; ++j)
    {
      elements[i].dof_values[j] = dof_values[global_dof_index];
      global_dof_index++;
    }

    global_dof_index--; //elements share a common boundary
  }
    



}


double PiecewisePolynomial::getValue(const double x)
{
  //check validity range
  if (x > dof_vertices.front() || x < dof_vertices.back())
  {
    std::cout << "Requested x value outside of domain of the polynomial!\n";
    
    return 0.0;
  }
  
  //first, we check if x is a global DOF
  for (size_t i=0; i<nb_dof; ++i)
    if (dof_vertices[i] == x) return dof_values[i];


  unsigned int element_index = 0;
  
  //if not, find the element it is in
  for (size_t i=0; i<nb_elements; ++i)
    if (elements[i].dof_vertices.front() > x && elements[i].dof_vertices.back() < x)
    {
      element_index = i;
      break;
    }

  //get the value from the corresponding element
  return elements[element_index].getValue(x);
}





Element::Element(const double edges [2], const size_t order)
{
  
  reference_vertices = quadrature_nodes[order];

  nb_dof = reference_vertices.size();

  dof_values.assign(nb_dof, 0.0);
  dof_vertices.resize(nb_dof);

  for (size_t i=0; i<nb_dof; ++i)
    dof_vertices[i] = referenceElementMap(reference_vertices[i], edges[0], edges[1]);
}


//Lagrange base funcion on the reference element
double Element::lagrangeBase(const double r, const unsigned int i)
{
  
  double l = 1;

  for (size_t j=0; j<nb_dof; ++j)
  {
     if (i == j) continue;

     l *= (r - reference_vertices[j]) / (reference_vertices[i] - reference_vertices[j]);

  }


  return l;
}


double Element::getValue(const double x)
{

  //coordinate on the reference element
  double r = realElementMap(x, dof_vertices.front(), dof_vertices.back());

  double y = 0;

  for (size_t i=0; i<nb_dof; ++i)
    y += dof_values[i] * lagrangeBase(r, i);  
    
    

  return y;
}



}