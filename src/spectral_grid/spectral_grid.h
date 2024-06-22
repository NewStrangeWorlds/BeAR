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


#ifndef _spectral_grid_h
#define _spectral_grid_h

#include <vector>
#include <string>

#include "spectral_band_type.h"


namespace bear {

//forward declaration
class GlobalConfig;
class Observation;



class SpectralGrid{
  public:
    SpectralGrid (GlobalConfig* global_config);
    ~SpectralGrid();

    std::vector<double> wavenumber_list;      //wavenumber list used to calculate the high-res spectra
    std::vector<double> wavelength_list;      //wavelength list used to calculate the high-res spectra

    double * wavenumber_list_gpu = nullptr;   //corresponding pointer to data on the GPU 
    double * wavelength_list_gpu = nullptr;

    std::vector<double> wavelengthToWavenumber(
      const std::vector<double>& wavelengths);
    std::vector<double> wavenumberToWavelength(
      const std::vector<double>& wavenumbers);

    double wavelengthToWavenumber(const double wavelength)
     {return 1.0/wavelength * 1e4;}
    double wavenumberToWavelength(const double wavenumber)
     {return 1.0/wavenumber * 1e4;}

    size_t nbSpectralPointsFull() {
      return nb_spectral_points_full;}
    size_t nbSpectralPoints() {
      return nb_spectral_points;}

    std::vector<size_t> spectralIndexList() {
      return index_list;}

    void sampleSpectralGrid(std::vector<Observation>& observations);

    void findBinEdges(
      const std::vector< std::vector<double> >& wavenumber_edges,
      std::vector<std::vector<size_t>>& edge_indices);

    size_t findClosestIndex(
      const double search_value,
      std::vector<double>& data,
      std::vector<double>::iterator it_start);

    std::vector<double> interpolateToWavenumberGrid(
      const std::vector<double>& data_x,
      const std::vector<double>& data_y,
      const bool log_interpolation);
    std::vector<double> interpolateToWavelengthGrid(
      const std::vector<double>& data_x,
      const std::vector<double>& data_y,
      const bool log_interpolation);
    std::vector<double> interpolateToWavelengthGrid(
      const std::vector<double>& data_x,
      const std::vector<double>& data_y,
      const std::vector<double>& new_x,
      const bool log_interpolation);
  private:
    GlobalConfig* config;

    std::vector<double> wavenumber_list_full; //the full, global wavenumber list, the opacities have been calculated at
    std::vector<double> wavelength_list_full; //the full, global wavelength list, the opacities have been calculated at

    std::vector<size_t> index_list;
    std::vector<std::vector<double>> observation_wavelength_edges;

    size_t nb_spectral_points_full;           //number of points in the global wavenumber list
    size_t nb_spectral_points;                //number of points in the spectral grid

    void loadWavenumberList();
    void createHeliosWavenumberList();

    void createHighResGrid(
      const std::vector<std::vector<size_t>>& edge_indices,
      std::vector<Observation>& observations);

    void createHighResGridConstWavenumber(
      const std::vector<std::vector<size_t>>& edge_indices,
      std::vector<int>& included_points);
    void createHighResGridConstWavelength(
      const std::vector<std::vector<size_t>>& edge_indices,
      std::vector<int>& included_points);
    void createHighResGridConstResolution(
      const std::vector<std::vector<size_t>>& edge_indices,
      std::vector<int>& included_points);

    void initDeviceMemory();

    size_t findClosestIndexDesc(
      const double search_value,
      std::vector<double>& data,
      std::vector<double>::iterator it_start);

    size_t findClosestIndexAsc(
      const double search_value,
      std::vector<double>& data,
      std::vector<double>::iterator it_start);
};



}


#endif
