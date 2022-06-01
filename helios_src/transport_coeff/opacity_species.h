

#ifndef OPACITY_SPECIES_H
#define OPACITY_SPECIES_H


#include <vector>
#include <string>
#include <iostream>

#include "../chemistry/chem_species.h"
#include "sampled_data.h"


namespace helios{


class GlobalConfig;
class SpectralGrid;



class OpacitySpecies {
  public:
    OpacitySpecies(const unsigned int index, const std::string name, const std::string folder) 
                  : species_index(index), species_name(name), species_folder(folder) 
                  {}
    virtual ~OpacitySpecies() {}
   
    bool dataAvailable() {return cross_section_available;}
    virtual void calcTransportCoefficients(const double temperature, const double pressure, const std::vector<double>& number_densities,
                                           std::vector<double>& absorption_coeff, std::vector<double>& scattering_coeff);
    virtual void calcTransportCoefficientsGPU(const double temperature, const double pressure, const std::vector<double>& number_densities,
                                              const size_t nb_grid_points, const size_t grid_point,
                                              double* absorption_coeff_device, double* scattering_coeff_device);
    
    const size_t species_index = 0;
    const std::string species_name = "";
    const std::string species_folder = "";
  protected:
    GlobalConfig* config;
    SpectralGrid* spectral_grid;

    double species_mass = 0;
    
    size_t pressure_reference_species = _TOTAL;
    std::vector<size_t> cia_collision_partner;

    bool cross_section_available = false;

    std::vector<SampledData> sampled_cross_sections;
    std::vector<std::vector<SampledData*>> ordered_data_list;

    void init();
    void orderDataList();

    virtual bool calcContinuumAbsorption(const double temperature, const std::vector<double>& number_densities, std::vector<double>& absorption_coeff) {return false;};
    virtual bool calcRalyleighCrossSections(std::vector<double>& cross_sections) {return false;};
    virtual void calcContinuumAbsorptionGPU(const double temperature, const std::vector<double>& number_densities,
                                            const size_t nb_grid_points, const size_t grid_point,
                                            double* absorption_coeff_device) {};

    void readFileList(const std::string file_path);
    
    std::vector<SampledData*> findClosestDataPoints(const double sampling_pressure, const double sampling_temperature);
    void checkDataAvailability(std::vector<SampledData*>& data_points);

    void calcAbsorptionCrossSections(const double local_pressure, const double local_temperature, std::vector<double>& cross_sections);
    bool calcScatteringCrossSections(std::vector<double>& cross_sections);

    void calcAbsorptionCoefficientsGPU(const double pressure, const double temperature, const double number_density,
                                       const size_t nb_grid_points, const size_t grid_point,
                                       double* absorption_coeff_device, double* scattering_coeff_device);
    
    double generalRayleighCrossSection(double reference_density, double refractive_index, double king_correction_factor, double wavenumber);
};


}

#endif
