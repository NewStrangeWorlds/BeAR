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


#ifndef _spectral_grid_h
#define _spectral_grid_h

#include <vector>
#include <string>

#include "spectral_band_type.h"


namespace helios {

//forward declaration
class GlobalConfig;



class SpectralGrid{
  public:
    SpectralGrid (GlobalConfig* global_config);
    ~SpectralGrid();

    void convertWavenumbersToWavelengths(const std::vector<double>& wavenumbers, std::vector<double>& wavelengths);
    void convertWavelengthsToWavenumbers(const std::vector<double>& wavelengths, std::vector<double>& wavenumbers);
    std::vector<double> convertWavelengthsToWavenumbers(const std::vector<double>& wavelengths);
    std::vector<double> convertWavenumbersToWavelengths(const std::vector<double>& wavenumbers);

    size_t nbSpectralPointsFull() {return nb_spectral_points_full;}
    size_t nbSpectralPoints() {return nb_spectral_points;}

    std::vector<size_t> spectralIndexList() {return index_list;}
    std::vector<double> wavenumberList(const std::vector<size_t>& indices);
    
    void sampleWavelengths(const std::vector< std::vector<double> >& band_edges, const double resolution, std::vector<size_t>& band_indices, std::vector<size_t>& bin_edge_indices);
    
    std::vector<size_t> sampleWavelengths(const std::vector<double>& wavelengths, const size_t nb_points_bin, std::vector<size_t>& edge_indices);
    std::vector<size_t> globalToLocalIndex(const std::vector<size_t>& global_indices);

    void findBinEdges(const std::vector< std::vector<double> >& wavenumber_edges, std::vector<double>& bin_edges, std::vector<size_t>& edge_indices);

    std::vector<double> interpolateToWavenumberGrid(const std::vector<double>& data_x, const std::vector<double>& data_y, const bool log_interpolation);
    std::vector<double> interpolateToWavelengthGrid(const std::vector<double>& data_x, const std::vector<double>& data_y, const bool log_interpolation);

    std::vector<double> wavenumber_list;                                         //wavenumber list used to calculate the high-res spectra
    std::vector<double> wavelength_list;                                         //wavelength list used to calculate the high-res spectra

    double * wavenumber_list_gpu = nullptr;                                      //corresponding pointer to data on the GPU 
    double * wavelength_list_gpu = nullptr;
  private:
    GlobalConfig* config;

    std::vector<double> wavenumber_list_full;                                    //the full, global wavenumber list, the opacities have been calculated at
    std::vector<double> wavelength_list_full;                                    //the full, global wavelength list, the opacities have been calculated at

    std::vector<size_t> index_list;

    size_t nb_spectral_points_full;                                              //number of points in the global wavenumber list
    size_t nb_spectral_points;                                                   //number of points in the spectral grid

    void loadWavenumberList();
    void addSampledIndices(const std::vector<size_t>& new_index_list);
};



}


#endif
