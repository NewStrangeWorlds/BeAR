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


#include "transport_coeff_single_species.h"


#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <omp.h>
#include <sstream>


#include <assert.h>
#include <cmath>

#include "../additional/aux_functions.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../CUDA_kernels/cross_section_kernels.h"


#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>



namespace helios{



void CIACoefficients::init(const std::string file_name, const std::vector<double>& wavenumber_grid, const bool use_gpu)
{

  readGeneralCIAFile(file_name, wavenumber_grid);



  if (use_gpu == true)
  {
    cross_sections_device.assign(nb_temperatures, nullptr);

    for (size_t i=0; i<nb_temperatures; ++i)
      moveToDevice(cross_sections_device[i], cross_sections[i]);

  }

}



CIACoefficients::~CIACoefficients()
{

  //in case we used the GPU, the data on it will be removed
  if (cross_sections_device.size() > 0)
  {

    for (size_t i=0; i<nb_temperatures; ++i)
      deleteFromDevice(cross_sections_device[i]);

  }


}




void CIACoefficients::readGeneralCIAFile(const std::string file_name, const std::vector<double>& wavenumber_grid)
{
  std::fstream file;

  file.open(file_name.c_str(), std::ios::in);
  assert (!file.fail( ));


  //read general header first
  std::string dummy_string;
  double dummy;

  size_t nb_wavenumbers;

  file >> dummy_string >> dummy >> dummy >> nb_wavenumbers;
  std::getline(file, dummy_string);

  std::cout << "number of wavenumbers: " << nb_wavenumbers << std::endl;

  std::vector<double> wavenumbers(nb_wavenumbers, 0.0);

  for (unsigned int i=0; i<nb_wavenumbers; i++)
    file >> wavenumbers[i] >> dummy;


  file.close();


  std::vector< std::vector<double> > cross_sections_file;

  //now read the data
  std::string header;
  file.open(file_name.c_str(), std::ios::in);

  while(std::getline(file, header))
  {
    if (file.eof()) break;

    std::stringstream  header_stream(header);
    double temperature;

    header_stream >> dummy_string >> dummy >> dummy >> dummy >> temperature;

    temperatures.push_back(temperature);

    cross_sections_file.resize(cross_sections_file.size() + 1);
    cross_sections_file.back().assign(nb_wavenumbers, 0.0);


    for (unsigned int j=0; j<nb_wavenumbers; j++)
    {
      std::string line;
      std::getline(file, line);
      std::stringstream line_stream(line);

      double data;

      line_stream >> dummy >> data;

      data = fabs(data);
      if (data < 1e-300) data = 1e-300;

      cross_sections_file.back()[j] = data;
    }


  }


  file.close();

  nb_temperatures = temperatures.size();

  std::cout << "number of temperatures: " << nb_temperatures << std::endl;


  //interpolate the cross sections from the file to the wavenumber grid
  cross_sections.resize(temperatures.size());


  for (size_t i=0; i<temperatures.size(); ++i)
  {

    cross_sections[i] = interpolateToWavenumberGrid(wavenumbers, cross_sections_file[i], wavenumber_grid, true);

    //later, we want to interpolate in log, so for performance reasons, we here calculate the log10 once and use that one
    for (size_t j=0; j<cross_sections[i].size(); ++j)
      cross_sections[i][j] = std::log10(cross_sections[i][j]);

  }

}



void CIACoefficients::findTemperatureIndices(const double temperature, size_t& index1, size_t& index2)
{

  if (nb_temperatures == 1 || temperature < temperatures.front())
  {
    index1 = 0;
    index2 = 0;

    return;
  }

  if (temperature > temperatures.back())
  {

    index1 = nb_temperatures-1;
    index2 = index1;

    return;
  }


  //find the two interpolation indices
  index1 = 0; index2 = 0;

  for (unsigned int i=1; i<nb_temperatures; i++)
  {
    if (temperature > temperatures[i-1] && temperature < temperatures[i])
    {
      index1 = i-1;
      index2 = i;

      break;
    }
    else if (temperature == temperatures[i])
    {

      index1 = i;
      index2 = i;

      break;
    }

  }



}



std::vector<double> CIACoefficients::temperatureInterpolation(const double temperature)
{

  size_t index1, index2;

  findTemperatureIndices(temperature, index1, index2);


  if (index1 == index2)
    return cross_sections[index1];


  std::vector<double> result(cross_sections[0].size(), 0.0);

  for (unsigned int i=0; i<cross_sections[0].size(); i++)
    result[i] = aux::linearInterpolation(temperatures[index1], temperatures[index2],
                                         cross_sections[index1][i], cross_sections[index2][i],
                                         temperature);


  return result;
}



std::vector<double> CIACoefficients::interpolateToWavenumberGrid(const std::vector<double>& wavenumber_data, const std::vector<double>& data,
                                                                 const std::vector<double>& wavenumber_grid,
                                                                 const bool interpolate_log)
{
  unsigned int nb_data_points = wavenumber_data.size();


  gsl_interp *spline = gsl_interp_alloc (gsl_interp_linear, nb_data_points);

  double *xa = new double[nb_data_points];
  double *ya = new double[nb_data_points];


  for (unsigned int i=0; i<nb_data_points; i++)
  {
    xa[i] = wavenumber_data[i];

    if (interpolate_log)
      ya[i] = log10(data[i]);
    else
      ya[i] = data[i];
  }



  gsl_interp_init (spline, xa, ya, nb_data_points);

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();



  std::vector<double> result(wavenumber_grid.size(), 0.0);


  for (unsigned int i=0; i<wavenumber_grid.size(); ++i)
  {
    if (wavenumber_grid[i] > wavenumber_data.back())
    {
      result[i] = 1e-300;

      continue;
    }

    if (wavenumber_grid[i] < wavenumber_data.front())
    {
      result[i] = 1e-300;

      continue;
    }

    if (interpolate_log)
      result[i] = pow(10., gsl_interp_eval (spline, xa, ya, wavenumber_grid[i],acc));
    else
      result[i] = gsl_interp_eval (spline, xa, ya, wavenumber_grid[i],acc);
  }


  gsl_interp_accel_free (acc);
  gsl_interp_free (spline);


  delete[] xa;
  delete[] ya;


  return result;
}




std::vector<double> CIACoefficients::calcCIACoefficients(const double temperature)
{

  return temperatureInterpolation(temperature);

}



void CIACoefficients::calcCIACoefficientsGPU(const double temperature, const double number_densities,
                                             const size_t nb_grid_points, const size_t grid_point,
                                             double* absorption_coeff_device)
{

  size_t index1, index2;

  findTemperatureIndices(temperature, index1, index2);


  calcCIACoefficientsHost(cross_sections_device[index1], cross_sections_device[index2],
                          temperatures[index1], temperatures[index2],
                          temperature, number_densities,
                          cross_sections[0].size(), nb_grid_points, grid_point,
                          absorption_coeff_device);


}



}
