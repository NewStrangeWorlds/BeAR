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


#ifndef _observations_h
#define _observations_h

#include <vector>
#include <string>
#include <iostream>

#include "../spectral_grid/spectral_grid.h"
#include "../config/global_config.h"
#include "../spectral_grid/spectral_band.h"


namespace bear {


//definition of different modifiers for an observation
//e.g. constants shifts
namespace observation_modifiers{
  enum id {none, shift_const}; 
  const std::vector<std::string> description {"none", "shift_const"};
}


//the class that describes an observation and its representation in form of a theoretical spectrum
//the class SpectralBand contains the computed represenation of the observation in the observational bands
//it also contains the required data and methods to convert a computed high-res spectrum into a simulated observation
class Observation{
  public:
    Observation(GlobalConfig* config_, SpectralGrid* spectral_grid_)
      : spectral_bands(config_, spectral_grid_)
      , config(config_)
      , spectral_grid(spectral_grid_)
      {}
    ~Observation();
    void init (
      const std::string& file_name,
      const std::string spectrum_modifier_id);
    std::string observationName() {return observation_name;}
    size_t nbPoints() {return data.size();}

    SpectralBands spectral_bands;                                       //representation of the theoretical spectrum in the observational bands
    std::vector<double> wavelength_edges = {0, 0};                      //the wavelength boundaries this observation needs from the high-res spectrum

    bool ascending_wavelengths = true;                                  //if the wavelength points of the orginal spectrum are in ascending order

    std::vector<double> data;                                           //observational data
    std::vector<double> data_error;                                     //observational error
    std::vector<double> likelihood_weight;                              //weight for the likelihood computation

    double* data_gpu = nullptr;                                         //pointer to the corresponding data on the GPU
    double* data_error_gpu = nullptr;                                   //pointer to the corresponding error on the GPU
    double* likelihood_weight_gpu = nullptr;                            //pointer to the corresponding likelihood weight on the GPU
    
    std::vector<double> instrument_profile_fwhm;                        //instrument profile (if it doesn't exist, the vector will have a size of 0)
    std::vector<double> filter_response;                                //filter response function
    std::vector<double> filter_response_weight;
    std::string filter_detector_type = "";
    double filter_response_normalisation = 0.0;

    double* filter_response_gpu = nullptr;
    double* filter_response_weight_gpu = nullptr;

    double* spectrum_filter_dev = nullptr;
    double* spectrum_convolved_dev = nullptr;

    void printObservationDetails();
    void setFilterResponseFunction();

    void initDeviceMemory();

    std::vector<double> applyFilterResponseFunction(const std::vector<double>& spectrum);
    void applyFilterResponseGPU(double* spectrum);

    std::vector<double> processModelSpectrum(
      const std::vector<double> spectrum, 
      const bool is_flux);
    void processModelSpectrumGPU(
      double* spectrum,
      double* spectrum_obs,
      const bool is_flux);

    void addShiftToSpectrumGPU(
      double* spectrum_obs,
      const double spectrum_shift);
    void addShiftToSpectrum(
      std::vector<double>& spectrum_bands,
      const double spectrum_shift);
    
    unsigned int nb_modifier_param = 0;
  private:
    std::string observation_name = "";
    GlobalConfig* config = nullptr;
    SpectralGrid* spectral_grid = nullptr;

    std::string filter_response_file_path = "";
    std::vector<std::vector<double>> filter_response_file;              //filter response function read from the file

    void loadFile(const std::string& file_name);

    std::vector<std::vector<double>> readFilterResponseFunction(const std::string& file_path);
    double filterResponseNormalisation(
      const std::vector<double>& filter_wavelength, const std::vector<double>& filter_response);

    bool readPhotometryData(std::fstream& file);
    bool readSpectroscopyData(std::fstream& file);
    bool readBandSpectroscopyData(std::fstream& file);

    observation_modifiers::id spectrum_modifier = observation_modifiers::id::none;

    void setSpectrumModifier(const std::string modifier_id);
};


}


#endif
