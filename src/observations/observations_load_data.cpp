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

  BandType band_type = PHOTOMETRY;
  std::string band_type_input;


  file >> band_type_input;


  if (band_type_input.at(0) == 'P' || band_type_input.at(0) == 'p')
    band_type = PHOTOMETRY;
  else if (band_type_input.at(0) == 'S' || band_type_input.at(0) == 's') 
    band_type = SPECTROSCOPY;
  else if (band_type_input.at(0) == 'B' || band_type_input.at(0) == 'b') 
    band_type = BAND_SPECTROSCOPY;
  else
  {
    std::string error_message = "Unsupported observation type *" + band_type_input + "* in file: " + file_name + "\n";
    throw InvalidInput(std::string ("Observation::loadFile"), error_message);
  }
    

  if (band_type == PHOTOMETRY) readPhotometryData(file);

  if (band_type == SPECTROSCOPY) readSpectroscopyData(file);

  if (band_type == BAND_SPECTROSCOPY) readBandSpectroscopyData(file);


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

  if (filter_response_file_path != "None" && filter_response_file_path != "none" && filter_response_file_path != "")
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
      std::string error_message = "Expected at least four data columns in observational file, but only found " + std::to_string(nb_data_columns) + "\n";
      throw InvalidInput(std::string ("Observation::readBandSpectroscopyData"), error_message);
    }
  }
  
  std::istringstream input(line);

  double lower_wavelength, upper_wavelength, flux_single, error_single, retrieval_weight;

  input >> lower_wavelength >> upper_wavelength >> flux_single >> error_single >> retrieval_weight;
  
  std::vector< std::vector<double> > wavelengths(1, {lower_wavelength, upper_wavelength});
 

  flux.push_back(flux_single);
  flux_error.push_back(error_single);

  if (nb_data_columns == 5)
    likelihood_weight.push_back(retrieval_weight);
  else
    likelihood_weight.push_back(1.0);


  //wavelengths should be ordered in descending order because the opacity data is organized in ascending wavenumbers
  if (wavelengths[0][0] < wavelengths[0][1])
    std::reverse(wavelengths[0].begin(), wavelengths[0].end());

  
  std::vector<double> band_centre(1, 0);
  band_centre[0] = wavelengths[0][0] - (wavelengths[0][0] - wavelengths[0][1]) * 0.5; 

  wavelength_edges[0] = wavelengths[0][0];
  wavelength_edges[1] = wavelengths[0][1];


  spectral_bands.init(wavelength_edges, wavelengths, band_centre, PHOTOMETRY);

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
  
  if (filter_response_file_path != "None" && filter_response_file_path != "none" && filter_response_file_path != "")
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
        std::string error_message = "Expected at least four data columns in observational file, but only found " + std::to_string(nb_data_columns) + "\n";
        throw InvalidInput(std::string ("Observation::readBandSpectroscopyData"), error_message);
      }
    }


    std::istringstream input(line);

    left_edge = 0;
    right_edge = 0;
    error = 0.0;
    line_profile_fwhm = 0.0;

    input >> left_edge >> right_edge >> flux_single >> error >> line_profile_fwhm >> retrieval_weight;

    if (left_edge != 0.0 && error != 0.0)
    {
      bin_edges.push_back({left_edge, right_edge});
      flux.push_back(flux_single);
      flux_error.push_back(error);
    }


    if (line_profile_fwhm != 0.0 && nb_data_columns > 4)
      instrument_profile_fwhm.push_back(line_profile_fwhm);

    if (nb_data_columns == 5)
      likelihood_weight.push_back(retrieval_weight);
    else
      likelihood_weight.push_back(1.0);
  }


  //wavelengths should be ordered in descending order because the opacity data is organized in ascending wavenumbers
  if (bin_edges.size() > 1)
  {
    if (bin_edges.front()[0] < bin_edges.back()[0])
    {
      std::reverse(bin_edges.begin(), bin_edges.end());
      std::reverse(flux.begin(), flux.end());
      std::reverse(flux_error.begin(), flux_error.end());
      std::reverse(instrument_profile_fwhm.begin(), instrument_profile_fwhm.end());
      std::reverse(likelihood_weight.begin(), likelihood_weight.end());
    }
  }


  for (size_t i=0; i<bin_edges.size(); ++i)
    if (bin_edges[i][0] < bin_edges[i][1]) std::reverse(bin_edges[i].begin(), bin_edges[i].end());


  std::vector<double> band_centres(bin_edges.size(), 0.0);

  for (size_t i=0; i<band_centres.size(); ++i)
    band_centres[i] = bin_edges[i][0] - (bin_edges[i][0] - bin_edges[i][1]) * 0.5; 


  wavelength_edges[0] = bin_edges.front()[0];
  wavelength_edges[1] = bin_edges.back()[1];


  spectral_bands.init(wavelength_edges, bin_edges, band_centres, BAND_SPECTROSCOPY);

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
  
  if (filter_response_file_path != "None" && filter_response_file_path != "none" && filter_response_file_path != "")
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
        std::string error_message = "Expected at least three data columns in observational file, but only found " + std::to_string(nb_data_columns) + "\n";
        throw InvalidInput(std::string ("Observation::readSpectroscopyData"), error_message);
      }
    }


    std::istringstream input(line);

    wavelength_single = 0.0;
    line_profile_fwhm = 0.0;
    retrieval_weight = 0.0;

    input >> wavelength_single >> flux_single >> error_single >> line_profile_fwhm >> retrieval_weight;

    if (wavelength_single != 0.0 && error_single != 0.0)
    {
      wavelengths.push_back(wavelength_single);
      flux.push_back(flux_single);
      flux_error.push_back(error_single);
    }


    if (line_profile_fwhm != 0.0 && nb_data_columns > 3)
      instrument_profile_fwhm.push_back(line_profile_fwhm);

    if (nb_data_columns == 5)
      likelihood_weight.push_back(retrieval_weight);
    else
      likelihood_weight.push_back(1.0);
  }


  //wavelengths should be ordered in descending order because the opacity data is organized in ascending wavenumbers
  if (wavelengths.size() > 1)
  {
    if (wavelengths[0] < wavelengths[1])
    {
      std::reverse(wavelengths.begin(), wavelengths.end());
      std::reverse(flux.begin(), flux.end());
      std::reverse(flux_error.begin(), flux_error.end());
      std::reverse(instrument_profile_fwhm.begin(), instrument_profile_fwhm.end());
      std::reverse(likelihood_weight.begin(), likelihood_weight.end());
    }
  }


  std::vector< std::vector<double> > bin_edges(wavelengths.size(), std::vector<double>(2, 0));


  //set up the bin edges
  for (size_t i=0; i<wavelengths.size(); ++i)
  {
    
    if (i == 0)
    {
      //first bin
      bin_edges.front()[1] = wavelengths[0] - (wavelengths[i] - wavelengths[i+1]) * 0.5;
      bin_edges.front()[0] = wavelengths[0] + (wavelengths[i] - wavelengths[i+1]) * 0.5;
    }
    else if (i == wavelengths.size()-1)
    {
      //last bin
      bin_edges.back()[1] = wavelengths[i] - (wavelengths[i-1] - wavelengths[i]) * 0.5;
      bin_edges.back()[0] = wavelengths[i] + (wavelengths[i-1] - wavelengths[i]) * 0.5;
    }
    else
    {
      bin_edges[i][0] = wavelengths[i-1] - (wavelengths[i-1] - wavelengths[i]) * 0.5;
      bin_edges[i][1] = wavelengths[i] - (wavelengths[i] - wavelengths[i+1]) * 0.5;
    }

  }

  //if we have an instrument profile, we extend the wavelength edges by 5 sigma
  if (instrument_profile_fwhm.size() != 0)
  {
    wavelength_edges[0] = bin_edges.front()[0] + 5.0 * instrument_profile_fwhm.front()/ 2.355;
    wavelength_edges[1] = bin_edges.back()[1] - 5.0 * instrument_profile_fwhm.back()/ 2.355;
  }
  else
  {
    wavelength_edges[0] = bin_edges.front()[0];
    wavelength_edges[1] = bin_edges.back()[1];
  }

  spectral_bands.init(wavelength_edges, bin_edges, wavelengths, SPECTROSCOPY);

  if (filter_response_file_path != "")
    filter_response_file = readFilterResponseFunction(filter_response_file_path);

  return true;
}



void Observation::printObservationDetails()
{
  std::cout << "\n";
  std::cout << "observation: " << observation_name << "\n";

  if (spectral_bands.bandType() == PHOTOMETRY)
  {
    std::cout << std::setprecision(5) << std::scientific << "observation type: photometry\n";
    std::cout << "Filter response function: " << filter_response_file_path << "\n";
    
    if (filter_response_file_path != "")
      std::cout << "Filter detector type: " << filter_detector_type << "\n";

    std::cout << "band edges: " << spectral_bands.edge_wavelengths[0][0] << "\t" << spectral_bands.edge_wavelengths[0][1] 
              << "\nband center: " << spectral_bands.center_wavelengths.front() << "\n";
    std::cout << "photometry flux: " << flux.front() << "\t error: " << flux_error.front() << "\n";
    std::cout << "likelihood weight: " << likelihood_weight.front() << "\n\n";
  }


  if (spectral_bands.bandType() == SPECTROSCOPY || spectral_bands.bandType() == BAND_SPECTROSCOPY)
  {
    if (spectral_bands.bandType() == BAND_SPECTROSCOPY)
      std::cout << "observation type: band-spectroscopy\n";
    else
     std::cout << "observation type: spectroscopy\n";
    
    if (filter_response_file_path != "")
    {
      std::cout << "Filter response function: " << filter_response_file_path << "\n";
      std::cout << "Filter detector type: " << filter_detector_type << "\n";
    }

    for (size_t i=0; i<flux.size(); ++i)
    {
      std::cout << std::setprecision(5) << std::scientific 
               << spectral_bands.center_wavelengths[i] << "\t" 
               << spectral_bands.edge_wavelengths[i][0] << "\t" 
               << spectral_bands.edge_wavelengths[i][1] 
               << "\t" << flux[i] << "\t" << flux_error[i];

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
