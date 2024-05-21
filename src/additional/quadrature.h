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


#ifndef _quadrature_h
#define _quadrature_h

#include <vector>
#include <stdexcept>

namespace helios{
namespace aux{



inline double quadratureTrapezoidal(const std::vector<double> &x, const std::vector<double> &y)
{

  if (x.size() != y.size())
    throw std::logic_error("trapezoidal quadrature: x and y must be the same size!\n");

  if (x.size() == 1)
    return y[0];


  double sum = 0.0;

  for (size_t i = 1; i < x.size(); ++i)
      sum += (x[i] - x[i-1]) * (y[i] + y[i-1]);


  return sum * 0.5;
}


inline double quadratureTrapezoidal(
  const std::vector<double> &x, 
  const std::vector<double> &y,
  const size_t idx_start,
  const size_t idx_end)
{

  if (x.size() != y.size())
    throw std::logic_error("trapezoidal quadrature: x and y must be the same size!\n");

  if (x.size() == 1)
    return y[0];


  double sum = 0.0;

  for (size_t i = idx_start+1; i < idx_end+1; ++i)
      sum += (x[i] - x[i-1]) * (y[i] + y[i-1]);


  return sum * 0.5;
}



}
}



#endif

