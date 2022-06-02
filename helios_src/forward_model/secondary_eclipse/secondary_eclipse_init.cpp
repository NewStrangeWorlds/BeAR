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


#include "secondary_eclipse.h"


#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>


#include "../../additional/exceptions.h"
#include "../../retrieval/retrieval.h"

#include "../../chemistry/select_chemistry.h"
#include "../../radiative_transfer/select_radiative_transfer.h"

#include "../../temperature/piecewise_poly_temperature.h"
#include "../../CUDA_kernels/data_management_kernels.h"


namespace helios{


//initialise radiative transfer model
void SecondaryEclipseModel::initRadiativeTransfer(const SecondaryEclipseConfig& model_config)
{
  radiative_transfer = selectRadiativeTransfer(model_config.radiative_transfer_model, 
                                               model_config.radiative_transfer_parameters, 
                                               model_config.nb_grid_points, 
                                               retrieval->config, &retrieval->spectral_grid);

}


//select and initialise the chemistry models
void SecondaryEclipseModel::initChemistry(const SecondaryEclipseConfig& model_config)
{
  chemistry.assign(model_config.chemistry_model.size(), nullptr);

  for (size_t i=0; i<model_config.chemistry_model.size(); ++i)
    chemistry[i] = selectChemistryModule(model_config.chemistry_model[i], model_config.chemistry_parameters[i], retrieval->config, model_config.atmos_boundaries);

  nb_total_chemistry_param = 0;

  for (auto & i : chemistry)
    nb_total_chemistry_param += i->nbParameters(); 
}



//select and initialise the chemistry models
void SecondaryEclipseModel::initTemperature(const SecondaryEclipseConfig& model_config)
{

  PiecewisePolynomialTemperature* temp = new PiecewisePolynomialTemperature(model_config.nb_temperature_elements, model_config.temperature_poly_degree, model_config.atmos_boundaries);
  temperature_profile = temp; 
  
}




//select and initialise the chemistry models
void SecondaryEclipseModel::initStellarSpectrum(const SecondaryEclipseConfig& model_config)
{
  std::fstream file;
  
  std::string file_path = retrieval->config->retrieval_folder_path + model_config.stellar_spectrum_file;

  file.open(file_path.c_str(), std::ios::in);


  if (file.fail())
    throw ExceptionFileNotFound(std::string ("SecondaryEclipseModel::initStellarSpectrum"), file_path);
 
  
  std::cout << "Reading stellar spectrum file " << file_path << "\n";
  
  std::vector<double> wavelength;  wavelength.reserve(5000000);
  std::vector<double> spectrum;  spectrum.reserve(5000000);

  std::string line;

  while (std::getline(file, line))
  {
    std::stringstream line_stream(line);

    double wavelength_in;
    double spectrum_in;

    if (!(line_stream >> wavelength_in >> spectrum_in)) continue;

    wavelength.push_back(wavelength_in);  
    spectrum.push_back(spectrum_in);
  }

  file.close();

  wavelength.shrink_to_fit(); spectrum.shrink_to_fit();


  //convert from W m-2 mu-1 to W m-2 cm
  for (size_t i=0; i<spectrum.size(); ++i)
    spectrum[i] = spectrum[i]*wavelength[i]*wavelength[i]/10000.;
    
  stellar_spectrum = retrieval->spectral_grid.interpolateToWavelengthGrid(wavelength, spectrum, false);


  binStellarSpectrum();
}



void SecondaryEclipseModel::binStellarSpectrum()
{
  stellar_spectrum_bands.assign(retrieval->observation_data.size(), 0.0);

  std::vector<double>::iterator it = stellar_spectrum_bands.begin();
  

  for (size_t i=0; i<retrieval->nb_observations; ++i)
  {
    std::vector<double> observation_bands;

    if (retrieval->observations[i].filter_response.size() == 0)
      retrieval->observations[i].spectral_bands.bandIntegrateSpectrum(stellar_spectrum, observation_bands);
    else
    {
      std::vector<double> filter_spectrum = retrieval->observations[i].applyFilterResponseFunction(stellar_spectrum);
      retrieval->observations[i].spectral_bands.bandIntegrateSpectrum(filter_spectrum, observation_bands);
    }
    
    std::copy(observation_bands.begin(), observation_bands.end(), it);
    it += observation_bands.size();
  }


  if (retrieval->config->use_gpu)
    moveToDevice(stellar_spectrum_bands_gpu, stellar_spectrum_bands);

}



}

