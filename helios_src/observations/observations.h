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


#ifndef _observations_h
#define _observations_h

#include <vector>
#include <string>



#include "../spectral_grid/spectral_band.h"


namespace helios {

//forward declaration
class Retrieval;


//the class that describes an observation and its representation in form of a theoretical spectrum
//the class SpectralBand contains the computed represenation of the observation in the observational bands
//it also contains the required data and methods to convert a computed high-res spectrum into a simulated observation
class Observation{
  public:
    ~Observation();
    void init (Retrieval* retrieval_ptr, const std::string& file_name);  //initialisation method that will read the file with the observational data  
    std::string observationName() const {return observation_name;} 

    SpectralBands spectral_bands;                                       //representation of the theoretical spectrum in the observational bands

    std::vector<double> flux;                                           //observational data
    std::vector<double> flux_error;                                     //observational error
    std::vector<double> instrument_profile_fwhm;                        //instrument profile (if it doesn't exist, the vector will have a size of 0)
    std::vector<double> filter_response;                                //filter response function
    std::vector<double> filter_response_weight;
    std::vector<double> likelihood_weight;                              //weight for the likelihood computation
    std::string filter_detector_type = "";
    double filter_response_normalisation = 0.0;
    
    double* filter_response_gpu = nullptr;
    double* filter_response_weight_gpu = nullptr;

    void printObservationDetails();
    void setFilterResponseFunction();
    std::vector<double> applyFilterResponseFunction(const std::vector<double>& spectrum);
  private:
    std::string observation_name = "";
    Retrieval* retrieval;

    std::string filter_response_file_path = "";
    std::vector<std::vector<double>> filter_response_file;              //filter response function read from the file

    void loadFile(const std::string& file_name);
    std::vector<std::vector<double>> readFilterResponseFunction(const std::string& file_path);
    double filterResponseNormalisation(const std::vector<double>& filter_wavelength, const std::vector<double>& filter_response);
    bool readPhotometryData(std::fstream& file);
    bool readSpectroscopyData(std::fstream& file);
    bool readBandSpectroscopyData(std::fstream& file);
};


}


#endif

