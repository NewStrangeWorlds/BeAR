/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2020 Daniel Kitzmann
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


#ifndef _piecewise_poly_h
#define _piecewise_poly_h


#include <cstddef>
#include <vector>



namespace helios {


class Element{
  public:
    Element(const double edges [2], const size_t order);
    std::vector<double> dof_vertices;
    std::vector<double> dof_values;
    size_t nb_dof = 0;

    double getValue(const double x);
  private:
    std::vector<double> reference_vertices; //the grid points on the reference element

    //maps the coordinate value r on the reference element [-1, +1] to the real element [x_l, x_r]
    double referenceElementMap(const double r, const double x_l, const double x_r) {return x_l + (1.0 + r)/2.0 * (x_r - x_l);}
    
    //maps the coordinate value x on the real element [x_l, x_r] to the reference element [-1, +1]
    double realElementMap(const double x, const double x_l, const double x_r) {return 2.0 * (x - x_l) / (x_r - x_l) - 1.0;}
    
    //the Lagrange base function on the reference elements
    double lagrangeBase(const double r, const unsigned int i);


};


class PiecewisePolynomial{
  public:
    PiecewisePolynomial(const size_t element_number, const size_t polynomial_order, const double domain_boundaries [2]);
    
    void setDOFvalues(const std::vector<double>& values);
    std::vector<double> getDOFvalues() {return dof_values;}

    //get the value of the polynomial at position x
    double getValue(const double x);

    size_t nb_dof = 0;
    std::vector<double> dof_vertices;
  private:
    std::vector<Element> elements;

    std::vector<double> dof_values; 

    size_t nb_elements = 0;
    size_t nb_edges = 0;
    unsigned int order = 0;
    std::vector<double> element_edges;
    
    size_t createElementGrid(const double domain_boundaries [2]);

    //maps the coordinate value r on the reference element [-1, +1] to the real element [x_l, x_r]
    double referenceElementMap(const double x, const double r, const double x_l, const double x_r) {return x_l + (1.0 + r)/2.0 * (x_r - x_l);}
};


}


#endif 
