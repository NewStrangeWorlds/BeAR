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


#include "transport_coeff_single_species.h"


#include "../config/global_config.h"
#include "../additional/physical_const.h"
#include "../chemistry/chem_species.h"

#include "../spectral_grid/spectral_grid.h"
#include "../CUDA_kernels/cross_section_kernels.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../additional/exceptions.h"


#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <omp.h>
#include <sstream>
#include <chrono>
#include <thread>
#include <cmath>


namespace helios{


void TransportCoefficientsSingleSpecies::init(std::vector<size_t>& sampling_index_list)
{
  sampling_points = sampling_index_list;
  nb_sampling_points = sampling_points.size();

  std::string file_path = config->cross_section_file_path;

  readFileList(file_path);

  sampled_cross_sections.resize(cross_section_data.size());


  for (unsigned int i=0; i<cross_section_data.size(); i++)
    sampled_cross_sections[i].init(cross_section_data[i].temperature,cross_section_data[i].pressure, config->use_gpu);
}



void TransportCoefficientsSingleSpecies::readFileList(const std::string& file_path)
{
  std::fstream file;
  std::string filelist_name = file_path+species_name+"/filelist.dat";

  file.open(filelist_name.c_str(), std::ios::in);

  if (file.fail())
  {
    std::cout << "\nNo line absorption cross-sections for species " << species_name << " found! Ignoring from here on.\n\n";
    cross_section_available = false;

    return;
  }
  else
    cross_section_available = true;


  std::string line;


  while (std::getline(file, line))
  {
    double t,p;
    std::string name;
    CrossSectionFile cross_section_file;

    std::istringstream input(line);

    input >> p >> t >> name;

    cross_section_file.pressure = p;
    cross_section_file.temperature = t;
    cross_section_file.filename = file_path+species_name+"/"+name;
    cross_section_file.is_loaded = false;

    cross_section_data.push_back(cross_section_file);
    
    //we briefly check if the file actually exists
    std::fstream opacity_file(cross_section_file.filename.c_str(), std::ios::binary | std::ios::in);
    
    if (opacity_file.fail())
    {
      std::string error_message = "file " + cross_section_file.filename + " of " + species_name + " from its filelist.dat not found!\n";
      throw ExceptionInvalidInput(std::string ("TransportCoefficientsSingleSpecies::readFileList"), error_message);
    }
      
    opacity_file.close();
  }


  file.close();


  //find min/max values for temperature and pressure
  for (size_t i=0; i<cross_section_data.size(); ++i)
  {
    if (max_temperature < cross_section_data[i].temperature) max_temperature = cross_section_data[i].temperature;
    if (max_pressure < cross_section_data[i].pressure) max_pressure = cross_section_data[i].pressure;
  }


  min_pressure = max_pressure;
  min_temperature = max_temperature;


  for (size_t i=0; i<cross_section_data.size(); ++i)
  {
    if (min_temperature > cross_section_data[i].temperature) min_temperature = cross_section_data[i].temperature;
    if (min_pressure > cross_section_data[i].pressure) min_pressure = cross_section_data[i].pressure;
  }  


  //for (unsigned int i=0; i<cross_section_data.size(); i++)
    //std::cout << i << "\t" << cross_section_data[i].filename << "\t" << cross_section_data[i].temperature << "\t" << cross_section_data[i].pressure << "\n";
}



//finds the four data points in the two-dimensional p-T grid of tabulated cross sections that are the closest to the given pressure and temperature
//this can probably be made more efficient by pre-ordering the list
void TransportCoefficientsSingleSpecies::findClosestDataPoints(unsigned int& lower_index1, unsigned int& lower_index2,
                                                               unsigned int& higher_index1, unsigned int& higher_index2,
                                                               const double pressure, const double temperature)
{
  double higher_pressure = max_pressure;
  double lower_pressure = min_pressure;
  double higher_temperature = min_temperature;
  double lower_temperature = min_temperature;


  double local_temperature = temperature;
  double local_pressure = pressure;




  for (unsigned int i=0; i<sampled_cross_sections.size(); i++)
  {
    if ((sampled_cross_sections[i].getTemperature() == local_temperature) && (sampled_cross_sections[i].getPressure() == local_pressure))
    {
      lower_index1 = i;
      higher_index1 = i;
      lower_index2 = i;
      higher_index2 = i;

      return;
    }
  }


  for (unsigned int i=0; i<sampled_cross_sections.size(); i++)
  {

    if (sampled_cross_sections[i].getTemperature() == local_temperature && sampled_cross_sections[i].getPressure() == local_pressure)
    {
      lower_index1 = i;
      higher_index1 = i;
      lower_index2 = i;
      higher_index2 = i;

      return;
    }

    if (sampled_cross_sections[i].getTemperature() <= local_temperature && sampled_cross_sections[i].getTemperature() >= lower_temperature &&
        sampled_cross_sections[i].getPressure() <= local_pressure && sampled_cross_sections[i].getPressure() >= lower_pressure)
    {
      lower_index1 = i;
      lower_temperature = sampled_cross_sections[i].getTemperature();
      lower_pressure = sampled_cross_sections[i].getPressure();
    }


    if (sampled_cross_sections[i].getTemperature() <= local_temperature && sampled_cross_sections[i].getTemperature() >= higher_temperature &&
        sampled_cross_sections[i].getPressure() >= local_pressure && sampled_cross_sections[i].getPressure() <= higher_pressure)
    {
      higher_index1 = i;
      higher_temperature = sampled_cross_sections[i].getTemperature();
      higher_pressure = sampled_cross_sections[i].getPressure();
    }

  }


  lower_temperature = max_temperature;
  higher_temperature = max_temperature;


  for (unsigned int i=0; i<sampled_cross_sections.size(); i++)
  {

    if (sampled_cross_sections[i].getPressure() == lower_pressure  && sampled_cross_sections[i].getTemperature() <= lower_temperature && sampled_cross_sections[i].getTemperature() >= local_temperature)
    {
      lower_index2 = i;
      lower_temperature = sampled_cross_sections[i].getTemperature();
    }


    if (sampled_cross_sections[i].getPressure() == higher_pressure  && sampled_cross_sections[i].getTemperature() <= higher_temperature && sampled_cross_sections[i].getTemperature() >= local_temperature)
    {
      higher_index2 = i;
      higher_temperature = sampled_cross_sections[i].getTemperature();
    }

  }


}



//calculation of the absorption coefficients for a specific temperature and pressure
//the method interpolates within a two-dimensional, tabulated cross section grid
void TransportCoefficientsSingleSpecies::calcAbsorptionCrossSections(const double pressure, const double temperature, std::vector<double>& cross_sections)
{
  cross_sections.assign(sampling_points.size(), 0.0);


  double local_temperature = temperature;

  if (local_temperature < min_temperature) local_temperature = min_temperature;
  if (local_temperature > max_temperature) local_temperature = max_temperature;


  double local_pressure = pressure;

  if (local_pressure < min_pressure) local_pressure = min_pressure;
  if (local_pressure > max_pressure) local_pressure = max_pressure;
  

  unsigned int lower_index1 = 0; 
  unsigned int higher_index1 = 0; 
  unsigned int lower_index2 = 0;
  unsigned int higher_index2 = 0;

  findClosestDataPoints(lower_index1,lower_index2,higher_index1,higher_index2,local_pressure,local_temperature);


  //std::cout << "pressure " << local_pressure << " temperature " << local_temperature << " total cross section number " << sampled_cross_sections.size() << std::endl;


  //std::cout << "lower index 1 " << lower_index1 << "\t" << sampled_cross_sections[lower_index1].getPressure() << "\t" << sampled_cross_sections[lower_index1].getTemperature() << "\t" << local_pressure << "\t" << local_temperature << std::endl;
  if (sampled_cross_sections[lower_index1].getSamplingStatus() == false)
    sampled_cross_sections[lower_index1].sampleCrossSections(cross_section_data[lower_index1], sampling_points);

  //std::cout << "lower index 2 " << lower_index2 << "\t" << sampled_cross_sections[lower_index2].getPressure() << "\t" << sampled_cross_sections[lower_index2].getTemperature() << std::endl;
  if (sampled_cross_sections[lower_index2].getSamplingStatus() == false)
    sampled_cross_sections[lower_index2].sampleCrossSections(cross_section_data[lower_index2], sampling_points);

  //std::cout << "higher index 1 " << higher_index1 << "\t" << sampled_cross_sections[higher_index1].getPressure() << "\t" << sampled_cross_sections[higher_index1].getTemperature() << std::endl;
  if (sampled_cross_sections[higher_index1].getSamplingStatus() == false)
    sampled_cross_sections[higher_index1].sampleCrossSections(cross_section_data[higher_index1], sampling_points);

  //std::cout << "higher index 2 " << higher_index2 << "\t" << sampled_cross_sections[higher_index2].getPressure() << "\t" << sampled_cross_sections[higher_index2].getTemperature() << std::endl;
  if (sampled_cross_sections[higher_index2].getSamplingStatus() == false)
    sampled_cross_sections[higher_index2].sampleCrossSections(cross_section_data[higher_index2], sampling_points);


  cross_sections.resize(sampling_points.size());



  std::vector<double> cross_sections_lower_index1;
  std::vector<double> cross_sections_lower_index2;
  std::vector<double> cross_sections_higher_index1;
  std::vector<double> cross_sections_higher_index2;

  sampled_cross_sections[higher_index1].getCrossSections(cross_sections_higher_index1);
  sampled_cross_sections[higher_index2].getCrossSections(cross_sections_higher_index2);
  sampled_cross_sections[lower_index1].getCrossSections(cross_sections_lower_index1);
  sampled_cross_sections[lower_index2].getCrossSections(cross_sections_lower_index2);



  if (lower_index1 != lower_index2)
  #pragma omp parallel for
    for (unsigned int i=0; i<cross_sections_lower_index1.size(); i++)
      cross_sections_lower_index1[i] = pow(10,linearInterpolation(sampled_cross_sections[lower_index1].getTemperature(),
                                                                  sampled_cross_sections[lower_index2].getTemperature(),
                                                                  log10(cross_sections_lower_index1[i]),
                                                                  log10(cross_sections_lower_index2[i]),
                                                                  local_temperature) );

  if (higher_index1 != higher_index2)
  #pragma omp parallel for
    for (unsigned int i=0; i<cross_sections_higher_index1.size(); i++)
      cross_sections_higher_index1[i] = pow(10,linearInterpolation(sampled_cross_sections[higher_index1].getTemperature(),
                                                                   sampled_cross_sections[higher_index2].getTemperature(),
                                                                   log10(cross_sections_higher_index1[i]),
                                                                   log10(cross_sections_higher_index2[i]),
                                                                   local_temperature) );

  if (lower_index1 != higher_index1)
  #pragma omp parallel for
    for (unsigned int i=0; i<cross_sections_lower_index1.size(); i++)
      cross_sections_lower_index1[i] = pow(10,linearInterpolation(log10(sampled_cross_sections[lower_index1].getPressure()),
                                                                  log10(sampled_cross_sections[higher_index1].getPressure()),
                                                                  log10(cross_sections_lower_index1[i]),
                                                                  log10(cross_sections_higher_index1[i]),
                                                                  log10(local_pressure)) );


  cross_sections = cross_sections_lower_index1;
}



//calculation of the absorption coefficients for a specific temperature and pressure
//the method interpolates within a two-dimensional, tabulated cross section grid
//this is version for doing the calculation on the GPU
//the four data points for the interpolation are still obtained on the CPU and then passed to the GPU
void TransportCoefficientsSingleSpecies::calcAbsorptionCoefficientsGPU(const double pressure, const double temperature, const double number_density,
                                                                       const size_t nb_grid_points, const size_t grid_point,
                                                                       double* absorption_coeff_device, double* scattering_coeff_device)
{
  double local_temperature = temperature;

  if (local_temperature < min_temperature) local_temperature = min_temperature;
  if (local_temperature > max_temperature) local_temperature = max_temperature;


  double local_pressure = pressure;

  if (local_pressure < min_pressure) local_pressure = min_pressure;
  if (local_pressure > max_pressure) local_pressure = max_pressure;


  unsigned int lower_index1 = 0; 
  unsigned int higher_index1 = 0; 
  unsigned int lower_index2 = 0;
  unsigned int higher_index2 = 0;

  findClosestDataPoints(lower_index1,lower_index2,higher_index1,higher_index2,local_pressure,local_temperature);


  if (sampled_cross_sections[lower_index1].getSamplingStatus() == false)
    sampled_cross_sections[lower_index1].sampleCrossSections(cross_section_data[lower_index1], sampling_points);

  if (sampled_cross_sections[lower_index2].getSamplingStatus() == false)
    sampled_cross_sections[lower_index2].sampleCrossSections(cross_section_data[lower_index2], sampling_points);

  if (sampled_cross_sections[higher_index1].getSamplingStatus() == false)
    sampled_cross_sections[higher_index1].sampleCrossSections(cross_section_data[higher_index1], sampling_points);

  if (sampled_cross_sections[higher_index2].getSamplingStatus() == false)
    sampled_cross_sections[higher_index2].sampleCrossSections(cross_section_data[higher_index2], sampling_points);


  calcCrossSectionsHost(sampled_cross_sections[lower_index1].cross_sections_device,
                        sampled_cross_sections[lower_index2].cross_sections_device,
                        sampled_cross_sections[higher_index1].cross_sections_device,
                        sampled_cross_sections[higher_index2].cross_sections_device,
                        sampled_cross_sections[lower_index1].getTemperature(),
                        sampled_cross_sections[lower_index2].getTemperature(),
                        sampled_cross_sections[lower_index1].getPressure(),
                        sampled_cross_sections[higher_index1].getPressure(),
                        local_temperature, local_pressure, number_density,
                        nb_sampling_points, nb_grid_points, grid_point,
                        absorption_coeff_device, scattering_coeff_device);
}





double TransportCoefficientsSingleSpecies::linearInterpolation(const double x1, const double x2, const double y1, const double y2, const double x)
{

  return y1 + (y2 - y1) * (x - x1)/(x2 - x1);

}



void TransportCoefficientsSingleSpecies::prepareCalculation(const double temperature, const double pressure)
{

  if (cross_section_available == true)
    checkDataAvailability(pressure, temperature);

}



//checks if all sampled cross sections needed later are available and samples them if not...
//currently not used
void TransportCoefficientsSingleSpecies::checkDataAvailability(const double pressure, const double temperature)
{
  unsigned int lower_index1 = 0; 
  unsigned int higher_index1 = 0; 
  unsigned int lower_index2 = 0;
  unsigned int higher_index2 = 0;


  findClosestDataPoints(lower_index1,lower_index2,higher_index1,higher_index2,pressure,temperature);


  //#pragma omp parallel sections num_threads(4)
  {
      //#pragma omp section
      {
        if (sampled_cross_sections[lower_index1].getSamplingStatus() == false)
          sampled_cross_sections[lower_index1].sampleCrossSections(cross_section_data[lower_index1], sampling_points);
      }

      //#pragma omp section
      {
        if (sampled_cross_sections[lower_index2].getSamplingStatus() == false && lower_index1 != lower_index2)
          sampled_cross_sections[lower_index2].sampleCrossSections(cross_section_data[lower_index2], sampling_points);
      }

      //#pragma omp section
      {
        if (sampled_cross_sections[higher_index1].getSamplingStatus() == false && lower_index1 != higher_index1 && lower_index2 != higher_index1)
          sampled_cross_sections[higher_index1].sampleCrossSections(cross_section_data[higher_index1], sampling_points);
      }


      //#pragma omp section
      {
        if (sampled_cross_sections[higher_index2].getSamplingStatus() == false && lower_index1 != higher_index2 && lower_index2 != higher_index2 && higher_index1 != higher_index2)
          sampled_cross_sections[higher_index2].sampleCrossSections(cross_section_data[higher_index2], sampling_points);
      }
  }

}




void TransportCoefficientsSingleSpecies::calcTransportCoefficients(const double temperature, const double pressure, const std::vector<double>& number_densities,
                                                                   std::vector<double>& absorption_coeff, std::vector<double>& scattering_coeff)
{
  if (number_densities[species_index] == 0) return;

  double gas_pressure = pressure;

  //K and Na are tabulated with respect to the H2 partial pressure
  if (species_name == "K" || species_name == "Na") gas_pressure = number_densities[_H2] * constants::boltzmann_k * temperature / 1e6;


  std::vector<double> cross_sections(nb_sampling_points, 0);


  if (cross_section_available == true) calcAbsorptionCrossSections(gas_pressure, temperature, cross_sections);

  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t i=0; i<nb_sampling_points; ++i)
    absorption_coeff[i] += cross_sections[i] * number_densities[species_index];


  //cross_sections.assign(nb_sampling_points, 0);

  /*calcScatteringCrossSections(cross_sections);

  for (size_t i=0; i<nb_sampling_points; ++i)
    scattering_coeff[i] += cross_sections[i] * number_densities[species_index];*/


  cross_sections.assign(nb_sampling_points, 0);


  //std::cout << "continuum absorption" << std::endl;
  if (calcContinuumAbsorption(temperature, number_densities, cross_sections) == true)
    for (size_t i=0; i<nb_sampling_points; ++i)
      absorption_coeff[i] += cross_sections[i];
}



void TransportCoefficientsSingleSpecies::calcTransportCoefficientsGPU(const double temperature, const double pressure, const std::vector<double>& number_densities,
                                                                      const size_t nb_grid_points, const size_t grid_point,
                                                                      double* absorption_coeff_device, double* scattering_coeff_device)
{
  if (number_densities[species_index] == 0) return;
  
  double gas_pressure = pressure;

  //K and Na are tabulated with respect to the H2 partial pressure
  if (species_name == "K" || species_name == "Na") gas_pressure = number_densities[_H2] * constants::boltzmann_k * temperature / 1e6;
  
  if (cross_section_available == true) calcAbsorptionCoefficientsGPU(gas_pressure, temperature, number_densities[species_index],
                                                                     nb_grid_points, grid_point,
                                                                     absorption_coeff_device, scattering_coeff_device);


  calcContinuumAbsorptionGPU(temperature, number_densities, nb_grid_points, grid_point, absorption_coeff_device);
}






void TransportCoefficientsSingleSpecies::calcScatteringCrossSections(std::vector<double>& cross_sections)
{

  calcRalyleighCrossSections(cross_sections);

}




double TransportCoefficientsSingleSpecies::generalRayleighCrossSection(const double reference_density, const double refractive_index, const double king_correction_factor,
                                                                       const double wavenumber)
{
  
  return 24. * pow(constants::pi, 3) * pow(wavenumber, 4) / (reference_density*reference_density)
         * pow((refractive_index*refractive_index - 1) / (refractive_index*refractive_index + 2),2) * king_correction_factor;

}




void CrossSectionFile::loadFile()
{
  std::fstream file;


  file.open(filename.c_str(), std::ios::binary | std::ios::in);


  if (file.fail()) std::cout << "cross section file " << filename << " not found! This shouldn't actually happen :-/ \n";
  if (file.fail()) exit(-1);


  unsigned int nb_data_points;

  file.read((char *) &nb_data_points, sizeof(int));


  cross_sections.resize(nb_data_points);


  for (unsigned int i=0; i<nb_data_points; ++i)
  {
    float x;

    file.read((char *) &x, sizeof x);

    cross_sections[i] = x;
  }


  file.close();

  is_loaded = true;
}


void CrossSectionFile::unloadData()
{
  std::vector<double>().swap (cross_sections);

  is_loaded = false;
}



}

