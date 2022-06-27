/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2022 Daniel Kitzmann
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

#ifndef SAMPLED_DATA_H
#define SAMPLED_DATA_H


#include <vector>
#include <string>
#include <cmath>


namespace helios{


class CrossSectionFile {
  public:
    CrossSectionFile(const std::string filename_in, const bool log_data)
        : filename(filename_in)
        , is_data_log(log_data) 
        {}
    void loadFile();
    void unloadData();
    
    const std::string filename = "";

    const bool is_data_log = false;
    bool is_loaded = false;

    std::vector<double> cross_sections;
};



class SampledData{
  public:
    SampledData(
      const double temperature_data, 
      const double pressure_data, 
      const std::string file_name, 
      const bool log_data, 
      const bool gpu_avail)
        : pressure(pressure_data)
        , log_pressure(std::log10(pressure_data))
        , temperature(temperature_data)
        , data_file(file_name, log_data)
        , use_gpu(gpu_avail) 
        {}
    ~SampledData();
    void sampleCrossSections(
      const std::vector<size_t>& sampling_list, const double species_mass);
    void deleteSampledData();

    const double pressure = 0.0;
    const double log_pressure = 0.0;
    const double temperature = 0.0;
    
    bool is_sampled = false;

    double* cross_sections_device = nullptr; //pointer to the cross section data on the GPU
    std::vector<double> cross_sections;
  private:
    CrossSectionFile data_file;
    const bool use_gpu = false;
};


}

#endif
