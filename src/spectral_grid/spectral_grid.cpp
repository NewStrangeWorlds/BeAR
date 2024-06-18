/*
* This file is part of the BeAR code (https://github.com/newstrangeworlds/BeAR).
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


#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "spectral_grid.h"

#include "spectral_band_type.h"
#include "../config/global_config.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../additional/exceptions.h"
#include "../observations/observations.h"
#include "spectral_band.h"


namespace helios{


SpectralGrid::SpectralGrid(GlobalConfig* global_config)
{
  config = global_config;

  loadWavenumberList();

  wavelength_list_full = wavenumberToWavelength(wavenumber_list_full);
}




void SpectralGrid::loadWavenumberList()
{
  std::string file_name = config->wavenumber_file_path;


  std::fstream file;

  file.open(file_name.c_str(), std::ios::in);


  if (file.fail())
    throw FileNotFound(std::string ("SpectralGrid::loadWavenumberList"), file_name);


  size_t nb_wavenumbers;
  file >> nb_wavenumbers;

  wavenumber_list_full.resize(nb_wavenumbers);


  for (std::vector<double>::iterator it = wavenumber_list_full.begin(); it != wavenumber_list_full.end(); ++it)
    file >> *it;


  file.close();
}



void SpectralGrid::sampleSpectralGrid(std::vector<Observation>& observations)
{
  std::vector<std::vector<double>> obs_wavenumber_edges;

  for (auto & o : observations)
    obs_wavenumber_edges.push_back(std::vector<double>{1.0/o.wavelength_edges[0]*1e4, 1.0/o.wavelength_edges[1]*1e4});

  std::sort(obs_wavenumber_edges.begin(), obs_wavenumber_edges.end(), [](std::vector<double> a, std::vector<double> b) {return a[0] < b[0];});

  std::vector<std::vector<double>> wavenumber_edges = 
    std::vector<std::vector<double>>{
      std::vector<double>{obs_wavenumber_edges[0][0], obs_wavenumber_edges[0][1]}
    };

  for (size_t i=1; i<obs_wavenumber_edges.size(); ++i)
  { 
    if (obs_wavenumber_edges[i][0] > wavenumber_edges.back()[0] && obs_wavenumber_edges[i][1] < wavenumber_edges.back()[1])
      continue;

    if (obs_wavenumber_edges[i][0] < wavenumber_edges.back()[1] && obs_wavenumber_edges[i][1] > wavenumber_edges.back()[1])
      wavenumber_edges.back()[1] = obs_wavenumber_edges[i][1];
    else
      wavenumber_edges.push_back(std::vector<double>{obs_wavenumber_edges[i][0], obs_wavenumber_edges[i][1]});
  }


  std::vector<std::vector<size_t>> edge_indices;
  findBinEdges(wavenumber_edges, edge_indices);


  createHighResGrid(edge_indices, observations);

  for (auto & o : observations)
    o.spectral_bands.setBandEdgeIndices(wavenumber_list);
}



void SpectralGrid::createHighResGrid(
  const std::vector<std::vector<size_t>>& edge_indices,
  std::vector<Observation>& observations)
{
  std::vector<int> included_points(wavenumber_list_full.size(), 0);

  for (size_t i=0; i<edge_indices.size(); ++i)
  {
    included_points[edge_indices[i][0]] = 1;

    size_t last_index = edge_indices[i][0];

    for (size_t j=edge_indices[i][0]; j<edge_indices[i][1]; ++j)
    {
      if (wavenumber_list_full[last_index] + config->spectral_resolution == wavenumber_list_full[j] 
        || (wavenumber_list_full[j] < wavenumber_list_full[last_index] + config->spectral_resolution && wavenumber_list_full[j+1] > wavenumber_list_full[last_index]+ config->spectral_resolution))
      { 
        included_points[j] = 1;
        last_index = j;
      }
      
      included_points[edge_indices[i][1]] = 1;
    }
  }

  //find the spectral indices of all band edges
  for (auto & o : observations)
  {
    auto it_start = wavenumber_list_full.begin();

    for (auto & b : o.spectral_bands.edge_wavenumbers)
    {
      const size_t idx_1 = findClosestIndex(b[0], wavenumber_list_full, it_start);
      const size_t idx_2 = findClosestIndex(b[1], wavenumber_list_full, it_start);

      included_points[idx_1] = 1;
      included_points[idx_2] = 1;
      
      it_start = wavenumber_list_full.begin() + idx_2;
    }
  }


  index_list.resize(0);
  index_list.reserve(wavelength_list_full.size());

  for (size_t i=0; i<included_points.size(); ++i)
  {
    if (included_points[i] == 1)
      index_list.push_back(i);
  }


  index_list.shrink_to_fit();

  wavenumber_list.assign(index_list.size(), 0);
  wavelength_list.assign(index_list.size(), 0);

  for (size_t i=0; i<index_list.size(); ++i)
  {
    wavenumber_list[i] = wavenumber_list_full[index_list[i]];
    wavelength_list[i] = wavelength_list_full[index_list[i]];
  }
  
  nb_spectral_points = wavelength_list.size();

  if (config->use_gpu)
    initDeviceMemory();
}


void SpectralGrid::initDeviceMemory()
{
  if (config->use_gpu == false)
    return;

  if (wavelength_list_gpu != nullptr)
  {
    deleteFromDevice(wavelength_list_gpu);
    deleteFromDevice(wavenumber_list_gpu);
  }

  moveToDevice(wavenumber_list_gpu, wavenumber_list);
  moveToDevice(wavelength_list_gpu, wavelength_list);
}


size_t SpectralGrid::findClosestIndexAsc(
  const double x,
  std::vector<double>& data,
  std::vector<double>::iterator start)
{
  auto iter_geq = std::lower_bound(
    start, 
    data.end(),
    x);

  if (iter_geq == start)
    return start - data.begin();

  double a = *(iter_geq - 1);
  double b = *(iter_geq);

  if (std::fabs(x - a) < fabs(x - b)) 
    return iter_geq - data.begin() - 1;

  return iter_geq - data.begin();
}


size_t SpectralGrid::findClosestIndexDesc(
  const double x,
  std::vector<double>& data,
  std::vector<double>::iterator start)
{
  auto iter_geq = std::lower_bound(
    start, 
    data.end(),
    x,
    std::greater<double>());

  if (iter_geq == start)
    return start - data.begin();

  double a = *(iter_geq - 1);
  double b = *(iter_geq);

  if (std::fabs(x - a) < fabs(x - b)) 
    return iter_geq - data.begin() - 1;

  return iter_geq - data.begin();
}


size_t SpectralGrid::findClosestIndex(
  const double x,
  std::vector<double>& data,
  std::vector<double>::iterator start)
{
  if (data.size() < 2)
    return 0;
  
  if (data[0] < data[1])
    return findClosestIndexAsc(x, data, start);

  if (data[0] > data[1])
    return findClosestIndexDesc(x, data, start);

  return 0;
}




void SpectralGrid::findBinEdges(
  const std::vector< std::vector<double> >& wavenumber_edges,
  std::vector<std::vector<size_t>>& edge_indices)
{
  edge_indices.assign(wavenumber_edges.size(), std::vector<size_t>{0, 0});

  for (size_t i=0; i<wavenumber_edges.size(); ++i)
  {
    auto it_left = std::lower_bound(
      wavenumber_list_full.begin(), 
      wavenumber_list_full.end(), 
      wavenumber_edges[i][0]);

    auto it_right = std::lower_bound(
      wavenumber_list_full.begin(), 
      wavenumber_list_full.end(), 
      wavenumber_edges[i][1]);

    edge_indices[i][0] = it_left - wavenumber_list_full.begin();

    if (wavenumber_list_full[edge_indices[i][0]] > wavenumber_edges[i][0])
      edge_indices[i][0] -= 1;

    edge_indices[i][1] = it_right - wavenumber_list_full.begin();
  }

}



/*std::vector<double> SpectralGrid::wavenumberList(const std::vector<size_t>& indices)
{
  std::vector<double> output(indices.size(), 0.0);

  for (size_t i=0; i<indices.size(); ++i)
    output[i] = wavenumber_list[indices[i]];

  return output;
}



std::vector<double> SpectralGrid::wavelengthList(const std::vector<size_t>& indices)
{
  std::vector<double> output(indices.size(), 0.0);

  for (size_t i=0; i<indices.size(); ++i)
    output[i] = wavelength_list[indices[i]];

  return output;
}*/



SpectralGrid::~SpectralGrid()
{
  if (config->use_gpu)
  {
    deleteFromDevice(wavelength_list_gpu);
    deleteFromDevice(wavenumber_list_gpu);
  }

}


}






