// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _aux_functions_h
#define _aux_functions_h

#include <vector>
#include <string>
#include "../chemistry/chem_species.h"

namespace bear{ namespace aux{


double planckFunctionWavenumber(const double temperature, const double wavenumber);
double linearInterpolation(const double x1, const double x2, const double y1, const double y2, const double x);
double triLinearInterpolation(const double x0, const double x1, 
                              const double y0, const double y1,
                              const double z0, const double z1, 
                              const double c000, const double c100, const double c010,
                              const double c110, const double c001, const double c101,
                              const double c011, const double c111,
                              const double x, const double y, const double z);
double voigtProfile(const double x, const double gaussian_width, const double lorentz_width);

double normalDistribution(const double mu, const double sigma, const double x);
double normalDistribution(const double sigma, const double x);

std::vector<double> interpolateToWavenumberGrid(const std::vector<double>& wavenumber_data, const std::vector<double>& data,
                                                const std::vector<double>& wavenumber_grid,
                                                const double outside_range_value,
                                                const bool interpolate_log);

std::vector<chemical_species_id> findChemicalSpecies(
  const std::vector<std::string>& species);

}}



#endif
