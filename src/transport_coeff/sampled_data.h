// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef SAMPLED_DATA_H
#define SAMPLED_DATA_H


#include <vector>
#include <string>
#include <cmath>


namespace bear{


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

    float* cross_sections_device = nullptr;
    std::vector<double> cross_sections;
  private:
    CrossSectionFile data_file;
    const bool use_gpu = false;
};


}

#endif
