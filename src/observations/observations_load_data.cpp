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


#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include "observations.h"

#include "../spectral_grid/spectral_band_type.h"
#include "../spectral_grid/spectral_band.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../additional/exceptions.h"


namespace bear{


void Observation::loadFile(const std::string& file_name)
{
  std::fstream file;

  file.open(file_name.c_str(), std::ios::in);


  if (file.fail())
    throw FileNotFound(std::string ("Observation::loadFile"), file_name);
 
  
  std::cout << "Reading observation file " << file_name << "\n";
  
  std::string line;
  std::getline(file, line);
  std::getline(file, line);
  
  std::getline(file, line);
  observation_name = line;
  
  std::getline(file, line);
  std::getline(file, line);

  std::string band_type_input;
  file >> band_type_input;
  
  band_type::id band_type = selectBandType(band_type_input, observation_name);


  if (band_type == band_type::id::photometry) readPhotometryData(file);

  if (band_type == band_type::id::spectroscopy) readSpectroscopyData(file);

  if (band_type == band_type::id::band_spectroscopy) readBandSpectroscopyData(file);


  file.close();
}




bool Observation::readPhotometryData(std::fstream& file)
{
  std::string line;

  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  
  //read the file name for the filter transmission function
  file >> filter_response_file_path;

  if (filter_response_file_path != "None" 
    && filter_response_file_path != "none" 
    && filter_response_file_path != "")
    filter_response_file_path = config->retrieval_folder_path + filter_response_file_path;
  else
    filter_response_file_path = "";


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);


  std::getline(file, line);
  
  //count the number of columns
  unsigned int nb_data_columns = 0;

  if (nb_data_columns == 0)
  {
    std::string item;

    std::stringstream ss(line);
    while (ss >> item ) nb_data_columns++;
      
    if (nb_data_columns < 4)
    {
      std::string error_message = 
        "Expected at least four data columns in observational file, but only found " 
        + std::to_string(nb_data_columns) + "\n";
      throw InvalidInput(
        std::string ("Observation::readBandSpectroscopyData"), 
        error_message);
    }
  }
  
  std::istringstream input(line);

  double lower_wavelength, upper_wavelength, flux_single, error_single, retrieval_weight;

  input >> lower_wavelength >> upper_wavelength >> flux_single >> error_single >> retrieval_weight;
  
  std::vector< std::vector<double> > wavelengths(1, {lower_wavelength, upper_wavelength});
 

  data.push_back(flux_single);
  data_error.push_back(error_single);

  if (nb_data_columns == 5)
    likelihood_weight.push_back(retrieval_weight);
  else
    likelihood_weight.push_back(1.0);

  ascending_wavelengths = areWavelengthsAscending(wavelengths);
  
  std::vector<double> band_centre = calcBinCenters(wavelengths);

  setObservationEdges(wavelengths);

  spectral_bands.init(
    wavelength_edges, 
    wavelengths, 
    band_centre, 
    band_type::id::photometry);

  if (filter_response_file_path != "")
    filter_response_file = readFilterResponseFunction(filter_response_file_path);

  return true;
}



bool Observation::readBandSpectroscopyData(std::fstream& file)
{
  double left_edge, right_edge, flux_single, error, line_profile_fwhm, retrieval_weight;
  std::vector< std::vector<double> > bin_edges; 

  std::string line;
  
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);

  //read the file name for the filter transmission function
  file >> filter_response_file_path;
  
  if (filter_response_file_path != "None" 
    && filter_response_file_path != "none" 
    && filter_response_file_path != "")
    filter_response_file_path = config->retrieval_folder_path + filter_response_file_path;
  else
    filter_response_file_path = "";

  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  
  unsigned int nb_data_columns = 0;

  while (std::getline(file, line))
  { 
    //count the number of columns
    if (nb_data_columns == 0)
    {
      std::string item;

      std::stringstream ss(line);
        while (ss >> item ) nb_data_columns++;
      
      if (nb_data_columns < 4)
      {
        std::string error_message = 
          "Expected at least four data columns in observational file, but only found " 
          + std::to_string(nb_data_columns) + "\n";
        throw InvalidInput(
          std::string ("Observation::readBandSpectroscopyData"), 
          error_message);
      }
    }


    std::istringstream input(line);

    left_edge = 0;
    right_edge = 0;
    error = 0.0;
    line_profile_fwhm = 0.0;

    input >> left_edge 
      >> right_edge 
      >> flux_single 
      >> error 
      >> line_profile_fwhm 
      >> retrieval_weight;

    if (left_edge != 0.0 && error != 0.0)
    {
      bin_edges.push_back({left_edge, right_edge});
      data.push_back(flux_single);
      data_error.push_back(error);
    }


    if (line_profile_fwhm != 0.0 && nb_data_columns > 4)
      instrument_profile_fwhm.push_back(line_profile_fwhm);

    if (nb_data_columns == 5)
      likelihood_weight.push_back(retrieval_weight);
    else
      likelihood_weight.push_back(1.0);
  }

  ascending_wavelengths = areWavelengthsAscending(bin_edges);

  std::vector<double> band_centres = calcBinCenters(bin_edges);

  setObservationEdges(bin_edges);

  spectral_bands.init(
    wavelength_edges, 
    bin_edges, 
    band_centres, 
    band_type::id::band_spectroscopy);

  if (filter_response_file_path != "")
    filter_response_file = readFilterResponseFunction(filter_response_file_path);

  return true;
}




bool Observation::readSpectroscopyData(std::fstream& file)
{
  std::vector<double> wavelengths;
  double wavelength_single, flux_single, error_single, line_profile_fwhm, retrieval_weight;
  std::string line;

  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  
  //read the file name for the filter transmission function
  file >> filter_response_file_path;
  
  if (filter_response_file_path != "None" 
    && filter_response_file_path != "none" 
    && filter_response_file_path != "")
    filter_response_file_path = config->retrieval_folder_path + filter_response_file_path;
  else
    filter_response_file_path = "";


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);

  unsigned int nb_data_columns = 0;

  while (std::getline(file, line))
  {
    //count the number of columns
    if (nb_data_columns == 0)
    {
      std::string item;

      std::stringstream ss(line);
        while (ss >> item ) nb_data_columns++;
      
      if (nb_data_columns < 3)
      {
        std::string error_message = 
          "Expected at least three data columns in observational file, but only found " 
          + std::to_string(nb_data_columns) + "\n";
        throw InvalidInput(
          std::string ("Observation::readSpectroscopyData"), 
          error_message);
      }
    }


    std::istringstream input(line);

    wavelength_single = 0.0;
    line_profile_fwhm = 0.0;
    retrieval_weight = 0.0;

    input >> wavelength_single 
      >> flux_single 
      >> error_single 
      >> line_profile_fwhm 
      >> retrieval_weight;

    if (wavelength_single != 0.0 && error_single != 0.0)
    {
      wavelengths.push_back(wavelength_single);
      data.push_back(flux_single);
      data_error.push_back(error_single);
    }


    if (line_profile_fwhm != 0.0 && nb_data_columns > 3)
      instrument_profile_fwhm.push_back(line_profile_fwhm);

    if (nb_data_columns == 5)
      likelihood_weight.push_back(retrieval_weight);
    else
      likelihood_weight.push_back(1.0);
  }
  
  ascending_wavelengths = areWavelengthsAscending(wavelengths);
  
  std::vector<std::vector<double>> bin_edges = calcBinEdges(wavelengths);
  
  setObservationEdges(bin_edges);
  
  spectral_bands.init(
    wavelength_edges, 
    bin_edges, 
    wavelengths, 
    band_type::id::spectroscopy);
  

  if (filter_response_file_path != "")
    filter_response_file = readFilterResponseFunction(filter_response_file_path);
  
  return true;
}



void Observation::printObservationDetails()
{
  std::cout << "\n";
  std::cout << "observation: " << observation_name << "\n";

  if (spectral_bands.bandType() == band_type::id::photometry)
  {
    std::cout << std::setprecision(5) << std::scientific << "observation type: photometry\n";
    std::cout << "Filter response function: " << filter_response_file_path << "\n";
    
    if (filter_response_file_path != "")
      std::cout << "Filter detector type: " << filter_detector_type << "\n";

    std::cout << "band edges: " << spectral_bands.edge_wavelengths[0][0] << "\t" << spectral_bands.edge_wavelengths[0][1] 
              << "\nband center: " << spectral_bands.center_wavelengths.front() << "\n";
    std::cout << "photometry flux: " << data.front() << "\t error: " << data_error.front() << "\n";
    std::cout << "likelihood weight: " << likelihood_weight.front() << "\n\n";
  }


  if (spectral_bands.bandType() == band_type::id::spectroscopy 
   || spectral_bands.bandType() == band_type::id::band_spectroscopy)
  {
    if (spectral_bands.bandType() == band_type::id::band_spectroscopy)
      std::cout << "observation type: band-spectroscopy\n";
    else
     std::cout << "observation type: spectroscopy\n";
    
    if (filter_response_file_path != "")
    {
      std::cout << "Filter response function: " << filter_response_file_path << "\n";
      std::cout << "Filter detector type: " << filter_detector_type << "\n";
    }

    for (size_t i=0; i<data.size(); ++i)
    {
      std::cout << std::setprecision(5) << std::scientific 
               << spectral_bands.center_wavelengths[i] << "\t" 
               << spectral_bands.edge_wavelengths[i][0] << "\t" 
               << spectral_bands.edge_wavelengths[i][1] 
               << "\t" << data[i] << "\t" << data_error[i];

      if (instrument_profile_fwhm.size() > 0)
        std::cout << "\t" << instrument_profile_fwhm[i];

      std::cout << "\t" << likelihood_weight[i] << "\n";
    }

    if (instrument_profile_fwhm.size() == 0)
      std::cout << "\n" << "No Instrument profile found.\n";

    std::cout << "\n";
  } 

}




}
