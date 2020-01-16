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


#ifndef TRANSPORT_COEFF_SINGLE_SPECIES_H
#define TRANSPORT_COEFF_SINGLE_SPECIES_H


#include <vector>
#include <string>
#include <iostream>


namespace helios{


class GlobalConfig;
class SpectralGrid;



class CIACoefficients{
  public:
    unsigned int species_id_1 = 0;
    unsigned int species_id_2 = 0;

    unsigned int nb_temperatures = 0;


    std::vector< std::vector<double> > cross_sections;
    std::vector<double> temperatures;

    //pointers to the cross sections on the gpu
    std::vector< double* > cross_sections_device;

    ~CIACoefficients();
    void init(const std::string file_name, const std::vector<double>& wavenumber_grid, const bool use_gpu);
    void readGeneralCIAFile(const std::string file_name, const std::vector<double>& wavenumber_grid);
    void findTemperatureIndices(const double temperature, size_t& index1, size_t& index2);
    std::vector<double> temperatureInterpolation(const double temperature);
    std::vector<double> interpolateToWavenumberGrid(const std::vector<double>& wavenumber_data, const std::vector<double>& data,
                                                    const std::vector<double>& wavenumber_grid,
                                                    const bool interpolate_log);
    std::vector<double> calcCIACoefficients(const double temperature);
    void calcCIACoefficientsGPU(const double temperature, const double number_densities,
                                const size_t nb_grid_points, const size_t grid_point,
                                double* absorption_coeff_device);
};



class CrossSectionFile {
  public:
    std::string filename;
    bool is_loaded;

    double pressure;
    double temperature;

    std::vector<double> cross_sections;

    void loadFile();
    void unloadData();
  private:

};



class SampledCrossSectionData{
  public:
    ~SampledCrossSectionData();
    void init(double temperature_data, double pressure_data);
    bool getSamplingStatus() {return is_sampled;}
    double getCrossSection(unsigned int wavelength_index);
    void getCrossSections(std::vector<double>& data_vector);
    double getTemperature() {return temperature;}
    double getPressure() {return pressure;}
    void sampleCrossSections(CrossSectionFile& input_data, const std::vector<size_t>& sampling_list);
    void deleteSampledData();

    double* cross_sections_device = nullptr; //pointer to the cross section data on the GPU
    std::vector<double> cross_sections;
  private:
    bool is_sampled;
    double pressure;
    double temperature;


};



class TransportCoefficientsSingleSpecies {
  public:
    virtual ~TransportCoefficientsSingleSpecies() {}
    void init(std::vector<size_t>& sampling_index_list);
    unsigned getSpeciesIndex() {return species_index;}
    std::string getSpeciesName() {return species_name;}
    virtual void prepareCalculation(const double temperature, const double pressure);
    virtual void calcTransportCoefficients(const double temperature, const double pressure, const std::vector<double>& number_densities,
                                           std::vector<double>& absorption_coeff, std::vector<double>& scattering_coeff);
    virtual void calcTransportCoefficientsGPU(const double temperature, const double pressure, const std::vector<double>& number_densities,
                                              const size_t nb_grid_points, const size_t grid_point,
                                              double* absorption_coeff_device, double* scattering_coeff_device);

  protected:
    GlobalConfig* config;
    SpectralGrid* spectral_grid;

    unsigned int species_index;
    std::string species_name;

    bool cross_section_available;

    std::vector<size_t> sampling_points;
    size_t nb_sampling_points;

    double min_temperature = 0.0;
    double max_temperature = 0.0;

    double max_pressure = 0.0;
    double min_pressure = 0.0;

    std::vector<CrossSectionFile> cross_section_data;
    std::vector<SampledCrossSectionData> sampled_cross_sections;


    virtual bool calcContinuumAbsorption(const double temperature, const std::vector<double>& number_densities, std::vector<double>& absorption_coeff) = 0;
    virtual void calcContinuumAbsorptionGPU(const double temperature, const std::vector<double>& number_densities,
                                            const size_t nb_grid_points, const size_t grid_point,
                                            double* absorption_coeff_device) = 0;
    virtual void calcRalyleighCrossSections(std::vector<double>& cross_sections) = 0;

    void readFileList(const std::string file_path);
    void findClosestDataPoints(unsigned int& lower_index1, unsigned int& lower_index2,
                               unsigned int& higher_index1, unsigned int& higher_index2,
                               double local_pressure, double local_temperature);
    void calcAbsorptionCrossSections(const double local_pressure, const double local_temperature, std::vector<double>& cross_sections);
    void calcAbsorptionCoefficientsGPU(const double pressure, const double temperature, const double number_density,
                                                                       const size_t nb_grid_points, const size_t grid_point,
                                                                       double* absorption_coeff_device, double* scattering_coeff_device);

    void calcScatteringCrossSections(std::vector<double>& cross_sections);

    double linearInterpolation(double x1, double x2, double y1, double y2, double x);
    void checkDataAvailability(double local_pressure, double local_temperature);

    double generalRayleighCrossSection(double reference_density, double refractive_index, double king_correction_factor, double wavenumber);
};



}


#endif
