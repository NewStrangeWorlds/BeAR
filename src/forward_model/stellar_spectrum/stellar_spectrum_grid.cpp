/*
* This file is part of the BeAR code (https://github.com/exoclime/BeAR).
* Copyright (C) 2024 Daniel Kitzmann
*
* BeAR is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* BeAR is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* BeAR directory under <LICENSE>. If not, see
* <http://www.gnu.org/licenses/>.
*/


#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <iomanip>
#include <sstream>
#include <algorithm> 

#include "stellar_spectrum_grid.h"

#include "../../spectral_grid/spectral_grid.h"
#include "../../additional/physical_const.h"
#include "../../additional/aux_functions.h"
#include "../../additional/exceptions.h"


namespace helios{


StellarSpectrumGrid::StellarSpectrumGrid(
  const std::string folder_path,
  SpectralGrid* spectral_grid_)
  : spectral_grid(spectral_grid_)
{
  nb_parameters = 3;

  std::string file_name = "grid_parameters.dat";
  
  std::string folder = folder_path;
  if (folder.back() != '/')
    folder += "+";

  std::string file_path = folder + file_name;

  std::vector<std::string> file_list = loadParameterFile(file_path);

  if (file_list.size() != effective_temperature.size()*log_g.size()*metallicity.size())
  {
    std::string error_message = "Stellar spectrum grid is inconsistent. The number of files does not correspond to the number of given effective temperatures, log g, and metallicities!\n";
    throw InvalidInput(std::string (file_path), error_message);
  }

  file_path = folder + wavelength_grid_file_name;

  loadWavelengthGrid(file_path);
  parseFileList(file_list, folder);
}


std::vector<std::string> StellarSpectrumGrid::loadParameterFile(
  const std::string& file_path)
{
  std::fstream file;

  file.open(file_path.c_str(), std::ios::in);

  if (file.fail())
    throw FileNotFound(std::string ("StellarSpectrumGrid::loadParameterFile"), file_path);


  std::cout << "Reading stellar spectra grid parameter file " << file_path << "\n";

  effective_temperature.reserve(100);

  std::string line;
  std::getline(file, line);

  double temperature = 0;
  std::istringstream temperature_input(line);

  while (temperature_input >> temperature)
    effective_temperature.push_back(temperature);

  effective_temperature.shrink_to_fit();


  std::getline(file, line);

  double logg = 0;
  std::istringstream log_g_input(line);

  while (log_g_input >> logg)
    log_g.push_back(logg);

  log_g.shrink_to_fit();


  std::getline(file, line);

  double met = 0;
  std::istringstream metallicity_input(line);

  while (metallicity_input >> met)
    metallicity.push_back(met);

  metallicity.shrink_to_fit();

  file >> wavelength_grid_file_name;

  std::getline(file, line);
  
  std::vector<std::string> file_list;
  file_list.reserve(effective_temperature.size()*metallicity.size()*log_g.size());
  
  while (std::getline(file, line))
    file_list.push_back(line);

  return file_list;
}



void StellarSpectrumGrid::loadWavelengthGrid(const std::string& file_path)
{
  std::fstream file;

  file.open(file_path.c_str(), std::ios::binary | std::ios::in | std::ios::ate);

  if (file.fail())
    throw FileNotFound(std::string ("StellarSpectrumGrid::loadWavelengthGrid"), file_path);


  int nb_data_points = file.tellg()/sizeof(float);
  file.seekg(0, std::ios::beg);

  grid_wavelengths.resize(nb_data_points);

  for (int i=0; i<nb_data_points; ++i)
  {
    float x;

    file.read((char *) &x, sizeof x);

    grid_wavelengths[i] = x;
  }

  file.close();
}



void StellarSpectrumGrid::parseFileList(
  const std::vector<std::string>& file_list,
  const std::string& folder)
{
  grid_spectra.reserve(effective_temperature.size()*metallicity.size()*log_g.size());

  for (auto & f : file_list)
  {
    std::istringstream input(f);
    double teff, logg, met;
    std::string file_name;

    input >> teff >> logg >> met >> file_name;

    std::string file_path = folder + file_name;

    grid_spectra.push_back(SampledStellarSpectrum(teff, logg, met, file_path));
  }


  grid.resize(effective_temperature.size());

  for (auto & g_logg : grid)
  {
    g_logg.resize(log_g.size());

    for (auto & g_met : g_logg)
      g_met.assign(metallicity.size(), nullptr);
  }


  for (auto & s : grid_spectra)
  {
    auto it = std::find(effective_temperature.begin(), effective_temperature.end(), s.temperature); 

    if (it == effective_temperature.end())  
    {
      std::string error_message = "Temperature " + std::to_string(s.temperature) + " not found in stellar grid description!";
      throw InvalidInput(std::string("StellarSpectrumGrid::parseFileList"), error_message);
    }
    
    size_t t_index = it - effective_temperature.begin(); 


    it = std::find(metallicity.begin(), metallicity.end(), s.metallicity); 

    if (it == metallicity.end())  
    {
      std::string error_message = "Metallicity " + std::to_string(s.metallicity) + " not found in stellar grid description!";
      throw InvalidInput(std::string("StellarSpectrumGrid::parseFileList"), error_message);
    }
    
    size_t m_index = it - metallicity.begin(); 


    it = std::find(log_g.begin(), log_g.end(), s.log_g); 

    if (it == log_g.end())
    {
      std::string error_message = "Log g " + std::to_string(s.log_g) + " not found in stellar grid description!";
      throw InvalidInput(std::string("StellarSpectrumGrid::parseFileList"), error_message);
    }
    
    size_t g_index = it - log_g.begin(); 

    grid[t_index][g_index][m_index] = &s;
  }


  //check the grid for completeness
  for (size_t i=0; i<effective_temperature.size(); ++i)
    for (size_t j=0; j<log_g.size(); ++j)
      for (size_t k=0; k<metallicity.size(); ++k)
      {
        if (grid[i][j][k] == nullptr)
        {
          std::string error_message = "Stellar spectra grid seems to be incomplete. Possibly duplicates?";
          throw InvalidInput(std::string("StellarSpectrumGrid::parseFileList"), error_message);
        }
      }

}



std::pair<size_t, size_t> StellarSpectrumGrid::findClosestGridPoints(
  const std::vector<double>& grid_points, 
  const double point)
{
  const size_t nb_points = grid_points.size();
  
  //treat boundary cases
  if (point <= grid_points.front())
    return std::pair<size_t, size_t>{0, 0};

  if (point >= grid_points.back())
    return std::pair<size_t, size_t>{nb_points-1, nb_points-1};


  std::pair<size_t, size_t> indices = {0,0};

  auto const it = std::lower_bound(grid_points.begin(), grid_points.end(), point);
  indices.second = std::distance(grid_points.begin(), it);
  
  if (grid_points[indices.second] != point)
    indices.first = indices.second - 1;
  else
    indices.first = indices.second;

  return indices;
}




std::vector<double> StellarSpectrumGrid::calcFlux(
  const std::vector<double>& parameter)
{
  std::vector<double> flux(spectral_grid->nbSpectralPoints(), 0.0);

  double temperature_int = parameter[0];
  double log_g_int = parameter[1];
  double metallicity_int = parameter[2];

  if (temperature_int < effective_temperature.front()) temperature_int = effective_temperature.front();
  if (temperature_int > effective_temperature.back()) temperature_int = effective_temperature.back();

  if (log_g_int < log_g.front()) log_g_int = log_g.front();
  if (log_g_int > log_g.back()) log_g_int = log_g.back();

  if (metallicity_int < metallicity.front()) metallicity_int = metallicity.front();
  if (metallicity_int > metallicity.back()) metallicity_int = metallicity.back();


  std::pair<size_t, size_t> t_i = findClosestGridPoints(effective_temperature, temperature_int);
  std::pair<size_t, size_t> g_i = findClosestGridPoints(log_g, log_g_int);
  std::pair<size_t, size_t> m_i = findClosestGridPoints(metallicity, metallicity_int);


  grid[t_i.first][g_i.first][m_i.first]->sampleSpectrum(spectral_grid, grid_wavelengths, false);
  grid[t_i.second][g_i.first][m_i.first]->sampleSpectrum(spectral_grid, grid_wavelengths, false);
  grid[t_i.first][g_i.second][m_i.first]->sampleSpectrum(spectral_grid, grid_wavelengths, false);
  grid[t_i.second][g_i.second][m_i.first]->sampleSpectrum(spectral_grid, grid_wavelengths, false);
  grid[t_i.first][g_i.first][m_i.second]->sampleSpectrum(spectral_grid, grid_wavelengths, false);
  grid[t_i.second][g_i.first][m_i.second]->sampleSpectrum(spectral_grid, grid_wavelengths, false);
  grid[t_i.first][g_i.second][m_i.second]->sampleSpectrum(spectral_grid, grid_wavelengths, false);
  grid[t_i.second][g_i.second][m_i.second]->sampleSpectrum(spectral_grid, grid_wavelengths, false);


  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t i=0; i<spectral_grid->nbSpectralPoints(); ++i)
  { 
     flux[i] = aux::triLinearInterpolation(
      effective_temperature[t_i.first], effective_temperature[t_i.second],
      log_g[g_i.first], log_g[g_i.second],
      metallicity[m_i.first], metallicity[m_i.second],
      grid[t_i.first][g_i.first][m_i.first]->spectrum[i],    //c000
      grid[t_i.second][g_i.first][m_i.first]->spectrum[i],   //c100
      grid[t_i.first][g_i.second][m_i.first]->spectrum[i],   //c010
      grid[t_i.second][g_i.second][m_i.first]->spectrum[i],  //c110
      grid[t_i.first][g_i.first][m_i.second]->spectrum[i],   //c001
      grid[t_i.second][g_i.first][m_i.second]->spectrum[i],  //c101
      grid[t_i.first][g_i.second][m_i.second]->spectrum[i],  //c011
      grid[t_i.second][g_i.second][m_i.second]->spectrum[i], //c111
      temperature_int, log_g_int, metallicity_int);
  }

  return flux;
}


} 
