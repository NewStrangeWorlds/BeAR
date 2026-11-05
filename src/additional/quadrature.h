// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _quadrature_h
#define _quadrature_h

#include <vector>
#include <stdexcept>

namespace bear{
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

