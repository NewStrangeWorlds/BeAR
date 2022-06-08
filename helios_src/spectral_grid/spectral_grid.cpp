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


#include "spectral_grid.h"


#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>


#include "spectral_band_type.h"
#include "../config/global_config.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../additional/exceptions.h"


namespace helios{


SpectralGrid::SpectralGrid(GlobalConfig* global_config)
{
  config = global_config;

  loadWavenumberList();

  convertWavenumbersToWavelengths(wavenumber_list_full, wavelength_list_full);
}




void SpectralGrid::loadWavenumberList()
{
  std::string file_name = config->wavenumber_file_path;


  std::fstream file;

  file.open(file_name.c_str(), std::ios::in);


  if (file.fail())
    throw ExceptionFileNotFound(std::string ("SpectralGrid::loadWavenumberList"), file_name);


  size_t nb_wavenumbers;
  file >> nb_wavenumbers;

  wavenumber_list_full.resize(nb_wavenumbers);


  for (std::vector<double>::iterator it = wavenumber_list_full.begin(); it != wavenumber_list_full.end(); ++it)
    file >> *it;


  file.close();
}



//identifies the indices from the opacity data related to the observational wavelengths
//distributes wavenumber points according to a given (constant) resolution
void SpectralGrid::sampleWavelengths(const std::vector< std::vector<double> >& wavelength_edges, const double resolution,
                                     std::vector<size_t>& band_indices, std::vector<size_t>& bin_edge_indices)
{
  std::vector< std::vector<double> > wavenumber_edges(wavelength_edges.size(), std::vector<double>(2, 0));

  
  for (size_t i=0; i<wavelength_edges.size(); ++i)
    convertWavelengthsToWavenumbers(wavelength_edges[i], wavenumber_edges[i]);

  
  std::vector<double> bin_edges;
  std::vector<size_t> edge_indices;
  
  //identify bind edges
  findBinEdges(wavenumber_edges, bin_edges, edge_indices);

  size_t nb_total_points_estimate = static_cast<size_t>((bin_edges.back() - bin_edges.front()) / resolution);
  
  band_indices.resize(0);
  band_indices.reserve(nb_total_points_estimate);

  
  const size_t nb_bins = edge_indices.size() - 1;

  bin_edge_indices.assign(nb_bins + 1, 0);


  //distribute the points within each bin
  for (size_t i=0; i<nb_bins; ++i)
  {
    size_t nb_bin_points = (bin_edges[i+1] - bin_edges[i]) / resolution;
    
    std::vector<size_t> single_bin_indices;
    single_bin_indices.reserve(nb_bin_points);

    
    double step = (edge_indices[i+1] - edge_indices[i]) / (1.0 * nb_bin_points-1);

    //in case there are less points available in a bin than desired, we simply take all points
    if (edge_indices[i+1] - edge_indices[i] < nb_bin_points)
    {
      for (size_t j=edge_indices[i]; j<edge_indices[i+1]+1; ++j)
        single_bin_indices.push_back(j);
    }
    else
    {
      single_bin_indices.push_back(edge_indices[i]);

      for (size_t j=1; j<nb_bin_points-1; ++j)
      {
        size_t index = single_bin_indices[0] + size_t (j * step);

        single_bin_indices.push_back(index);
      }


      single_bin_indices.push_back(edge_indices[i+1]);
    }


    //add the points of each bin to the band's list and identify the local indices of the bin edges
    if (i==0)
    {
      band_indices.insert(band_indices.end(), single_bin_indices.begin(), single_bin_indices.end());
      bin_edge_indices[0] = 0;
    }
    else
      band_indices.insert(band_indices.end(), single_bin_indices.begin()+1, single_bin_indices.end());


    bin_edge_indices[i+1] = band_indices.size()-1;
  }
  
  
  //add band indices to the global list
  addSampledIndices(band_indices);
}



void SpectralGrid::findBinEdges(const std::vector< std::vector<double> >& wavenumber_edges, std::vector<double>& bin_edges, std::vector<size_t>& edge_indices)
{
  bin_edges.assign(wavenumber_edges.size()+1, 0.0);
  edge_indices.assign(wavenumber_edges.size()+1, 0);


  for (size_t i=0; i<wavenumber_edges.size(); ++i)
    bin_edges[i] = wavenumber_edges[i][0];

  bin_edges.back() = wavenumber_edges.back()[1];


  //find the associated indices from the global list
  size_t j_start = 0;

  for (size_t i=0; i<bin_edges.size(); ++i)
  {

    for (size_t j=j_start; j<wavenumber_list_full.size(); ++j)
    {

      if ( wavenumber_list_full[j] < bin_edges[i] && wavenumber_list_full[j+1] > bin_edges[i] )
      {

        if ( std::abs(wavenumber_list_full[j] - bin_edges[i])  < std::abs(wavenumber_list_full[j+1] - bin_edges[i]) )
          edge_indices[i] = j;
        else
          edge_indices[i] = j+1;


        j_start = j;

        break;
      }

      if ( wavenumber_list_full[j] == bin_edges[i])
      {
        edge_indices[i] = j;

        j_start = j;

        break;
      }

    }


  }


  //change the bin edges to correspond to wavenumbers from those of the tabulated data
  for (size_t i = 0; i<bin_edges.size(); ++i)
    bin_edges[i] = wavenumber_list_full[edge_indices[i]];
}



//identifies the indices from the opacity data related to the observational wavelengths
std::vector<size_t> SpectralGrid::sampleWavelengths(const std::vector<double>& wavelengths, const size_t nb_points_bin, std::vector<size_t>& edge_indices)
{
  std::vector<double> wavenumbers;
  convertWavelengthsToWavenumbers(wavelengths, wavenumbers);


  std::vector<size_t> indices(wavelengths.size()*2, 0);


  size_t j_start = 0;

  for (size_t i=0; i<wavelengths.size(); ++i)
  {


    for (size_t j=j_start; j<wavelength_list_full.size(); ++j)
    {

      if ( (wavelength_list_full[j] > wavelengths[i] && wavelength_list_full[j+1] < wavelengths[i]))
      {
        indices[2*i] = j;
        indices[2*i+1] = j+1;

        j_start = j;

        break;
      }
      if ( wavelength_list_full[j] == wavelengths[i])
      {
        indices[2*i] = j;
        indices[2*i+1] = j;

        j_start = j;

        break;
      }

    }


  }

  addSampledIndices(indices);


  return indices;
}



//adds a new list of indices to the old one, orders it, and removes duplicates
void SpectralGrid::addSampledIndices(const std::vector<size_t>& new_index_list)
{ 
  index_list.insert(index_list.end(), new_index_list.begin(), new_index_list.end());


  std::sort (index_list.begin(), index_list.end());


  //remove duplicate indices
  //erase is not very economic but since we only do it once at the start probably OK
  std::vector<size_t>::iterator it = index_list.begin()+1;

  while (it != index_list.end())
  {
  
    if (*it == *(it-1))
    {
      it = index_list.erase(it);

    }

    else
      ++it;
  }


  wavenumber_list.assign(index_list.size(), 0);
  wavelength_list.assign(index_list.size(), 0);


  for (size_t i=0; i<index_list.size(); ++i)
  {
    wavenumber_list[i] = wavenumber_list_full[index_list[i]];
    wavelength_list[i] = wavelength_list_full[index_list[i]];
  }


  nb_spectral_points = index_list.size();


  //move the lists to the GPU, if necessary
  if (config->use_gpu)
  { 
    if (wavelength_list_gpu != nullptr)
    {
      deleteFromDevice(wavelength_list_gpu);
      deleteFromDevice(wavenumber_list_gpu);
    }

    moveToDevice(wavenumber_list_gpu, wavenumber_list);
    moveToDevice(wavelength_list_gpu, wavelength_list);
  }

}



//finds the local spectral index in the sampled list based on the global index
std::vector<size_t> SpectralGrid::globalToLocalIndex(const std::vector<size_t>& global_indices)
{
  std::vector<size_t> local_indices(global_indices.size(), 0);


  for (size_t i=0; i<global_indices.size(); ++i)
    for (size_t j=0; j<index_list.size(); ++j)
      if (global_indices[i] == index_list[j])
      {
        local_indices[i] = j;

        break;
      }


  return local_indices;
}



std::vector<double> SpectralGrid::wavenumberList(const std::vector<size_t>& indices)
{
  std::vector<double> output(indices.size(), 0.0);


  for (size_t i=0; i<indices.size(); ++i)
    output[i] = wavenumber_list[indices[i]];


  return output;
}




SpectralGrid::~SpectralGrid()
{
  if (config->use_gpu)
  {
    deleteFromDevice(wavelength_list_gpu);
    deleteFromDevice(wavenumber_list_gpu);
  }

}



}






