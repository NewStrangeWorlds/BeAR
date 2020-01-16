/*
* This file is part of the Helios-r2 code (https://github.com/exoclime/Helios-r2).
* Copyright (C) 2019 Daniel Kitzmann
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


#include "observations.h"


#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include "../spectral_grid/spectral_band_type.h"
#include "../spectral_grid/spectral_band.h"
#include "../retrieval/retrieval.h"
#include "../CUDA_kernels/data_management_kernels.h"


namespace helios{


//initialisation method that reads the observational data and sets up the basic data structures
//note that this might be better moved into the constructor
//we chose a separate method here because it can return a boolean that returns a success/failure
//within a constructor this could be achived via an exception but I was too lazy to implement this for now
//(might change later)
bool Observation::init (Retrieval* retrieval_ptr, const std::string file_name)
{
  retrieval = retrieval_ptr;


  bool file_loaded = loadFile(file_name);


  return file_loaded;
}


bool Observation::loadFile(const std::string file_name)
{
  std::fstream file;

  file.open(file_name.c_str(), std::ios::in);


  if (file.fail())
  {
    std::cout << "Could not open observation file: " << file_name << "\n";
    return false;
  }
  else
  {
    std::cout << "Reading observation file " << file_name << "\n";
  }
  
  

  std::string line;
  std::getline(file, line);
  std::getline(file, line);

  file >> observation_name;

  std::getline(file, line);
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
    std::cout << "Unsupported observation type \"" << band_type_input << "\" in file: " << file_name << "\n"; exit(0);
    return false;
  }


  if (band_type == PHOTOMETRY) readPhotometryData(file);

  if (band_type == SPECTROSCOPY) readSpectroscopyData(file);

  if (band_type == BAND_SPECTROSCOPY) readBandSpectroscopyData(file);


  file.close();

  
  return true;
}




bool Observation::readPhotometryData(std::fstream& file)
{
  std::string line;

  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);


  std::string filter_file;

  file >> filter_file;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);


  double lower_wavelength, upper_wavelength, flux_single, error_single;

  file >> lower_wavelength >> upper_wavelength >> flux_single >> error_single;
  
  std::vector< std::vector<double> > wavelengths(1, {lower_wavelength, upper_wavelength});
 

  flux.push_back(flux_single);
  flux_error.push_back(error_single);


  //wavelengths should be ordered in descending order because the opacity data is organized in ascending wavenumbers
  if (wavelengths[0][0] < wavelengths[0][1])
    std::reverse(wavelengths[0].begin(), wavelengths[0].end());

  
  std::vector<double> band_centre(1, 0);
  band_centre[0] = wavelengths[0][0] - (wavelengths[0][0] - wavelengths[0][1]) * 0.5; 

 
  spectral_bands.init(retrieval->config, &retrieval->spectral_grid, wavelengths, band_centre, PHOTOMETRY, filter_file);

  
  return true;
}



bool Observation::readBandSpectroscopyData(std::fstream& file)
{

  double left_edge, right_edge, flux_single, error, line_profile_fwhm;
  std::vector< std::vector<double> > bin_edges; 

    
  std::string line;

  while (std::getline(file, line))
  {
    std::istringstream input(line);

    left_edge = 0;
    right_edge = 0;
    error = 0.0;
    line_profile_fwhm = 0.0;

    input >> left_edge >> right_edge >> flux_single >> error >> line_profile_fwhm;

    if (left_edge != 0.0 && error != 0.0)
    {
      bin_edges.push_back({left_edge, right_edge});
      flux.push_back(flux_single);
      flux_error.push_back(error);    
    }


    if (line_profile_fwhm != 0.0)     
      instrument_profile_fwhm.push_back(line_profile_fwhm);
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
    }
  }


  for (size_t i=0; i<bin_edges.size(); ++i)
    if (bin_edges[i][0] < bin_edges[i][1]) std::reverse(bin_edges[i].begin(), bin_edges[i].end());


  std::string filter_file = "";



  std::vector<double> band_centres(bin_edges.size(), 0.0);

  for (size_t i=0; i<band_centres.size(); ++i)
    band_centres[i] = bin_edges[i][0] - (bin_edges[i][0] - bin_edges[i][1]) * 0.5; 

    
  spectral_bands.init(retrieval->config, &retrieval->spectral_grid, bin_edges, band_centres, SPECTROSCOPY, filter_file);


  return true;
}




bool Observation::readSpectroscopyData(std::fstream& file)
{

  std::vector<double> wavelengths, error_m;
  double wavelength_single, flux_single, error_p_single, error_m_single, line_profile_fwhm;

    
  std::string line;

  while (std::getline(file, line))
  {
    std::istringstream input(line);

    wavelength_single = 0.0;
    error_m_single = 0.0;
    line_profile_fwhm = 0.0;

    input >> wavelength_single >> flux_single >> error_p_single >> error_m_single >> line_profile_fwhm;

    if (wavelength_single != 0.0 && error_m_single != 0.0)
    {
      wavelengths.push_back(wavelength_single);
      flux.push_back(flux_single);
      flux_error.push_back(error_p_single);
      error_m.push_back(error_m_single);
    }


    if (line_profile_fwhm != 0.0)     
      instrument_profile_fwhm.push_back(line_profile_fwhm);
  }


  //wavelengths should be ordered in descending order because the opacity data is organized in ascending wavenumbers
  if (wavelengths.size() > 1)
  {
    if (wavelengths[0] < wavelengths[1])
    {
      std::reverse(wavelengths.begin(), wavelengths.end());
      std::reverse(flux.begin(), flux.end());
      std::reverse(flux_error.begin(), flux_error.end());
      std::reverse(error_m.begin(), error_m.end());
      std::reverse(instrument_profile_fwhm.begin(), instrument_profile_fwhm.end());
    }
  }


  std::string filter_file = "";



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
    

    
    
  spectral_bands.init(retrieval->config, &retrieval->spectral_grid, bin_edges, wavelengths, SPECTROSCOPY, filter_file);


  return true;
}



void Observation::printObservationDetails()
{

  std::cout << "\n";
  std::cout << "observation: " << observation_name << "\n";

  if (spectral_bands.bandType() == PHOTOMETRY)
  {
    std::cout << std::setprecision(5) << std::scientific << "observation type: photometry\n";
    std::cout << "band edges: " << spectral_bands.band_edges_wavelength[0][0] << "\t" << spectral_bands.band_edges_wavelength[0][1] 
              << "\nband center: " << spectral_bands.band_centers_wavelength.front() << "\n";
    std::cout << "photometry flux: " << flux.front() << "\t error: " << flux_error.front() << "\n\n";
  }


  if (spectral_bands.bandType() == SPECTROSCOPY || spectral_bands.bandType() == BAND_SPECTROSCOPY)
  {
    std::cout << "observation type: spectroscopy\n";

    for (size_t i=0; i<flux.size(); ++i)
    {
      std::cout << std::setprecision(5) << std::scientific 
               << spectral_bands.band_centers_wavelength[i] << "\t" 
               << spectral_bands.band_edges_wavelength[i][0] << "\t" 
               << spectral_bands.band_edges_wavelength[i][1] 
               << "\t" << flux[i] << "\t" << flux_error[i];

      if (instrument_profile_fwhm.size() > 0)
        std::cout << "\t" << instrument_profile_fwhm[i] << "\n";
      else
        std::cout << "\n";
    }

    if (instrument_profile_fwhm.size() == 0)
      std::cout << "\n" << "No Instrument profile found.\n";

    std::cout << "\n";
  } 
    


}




}
