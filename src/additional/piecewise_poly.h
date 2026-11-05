// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _piecewise_poly_h
#define _piecewise_poly_h


#include <cstddef>
#include <vector>



namespace bear {


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
    PiecewisePolynomial(
      const size_t element_number, 
      const size_t polynomial_order, 
      const std::vector<double>& domain_boundaries);
    
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
