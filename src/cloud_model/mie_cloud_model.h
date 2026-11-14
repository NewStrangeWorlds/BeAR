// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _mie_cloud_model_h
#define _mie_cloud_model_h


#include <vector>
#include <string>

#include "cloud_model.h"
#include "grey_cloud_model.h"
#include "../forward_model/atmosphere/atmosphere.h"


namespace bear {


//Non-grey cloud layer with optical properties taken from a table of Mie
//calculations for log-normal particle size distributions.
//The table provides the mean extinction cross section per particle, the
//single scattering albedo, and the asymmetry parameter as a function of
//the median particle radius r_g, the geometric standard deviation sigma_g,
//and wavelength. The cloud's optical depth at a reference wavelength is a
//free parameter; the wavelength dependence follows the extinction cross section.
class MieCloudModel: public GreyCloudModel{
  public:
    MieCloudModel(
      const std::vector<std::string>& parameters,
      const std::string& retrieval_folder);
    virtual ~MieCloudModel();
    virtual void opticalProperties(
      const std::vector<double>& parameters, 
      const Atmosphere& atmosphere,
      SpectralGrid* spectral_grid,
      std::vector<std::vector<double>>& optical_depth, 
      std::vector<std::vector<double>>& single_scattering, 
      std::vector<std::vector<double>>& asym_param);
    virtual void opticalPropertiesGPU(
      const std::vector<double>& parameters, 
      const Atmosphere& atmosphere,
      SpectralGrid* spectral_grid,
      float* optical_depth_dev, 
      float* single_scattering_dev, 
      float* asym_param);

    //interpolates the table to a given particle size and distribution width;
    //returns the log of the extinction cross section, the single scattering
    //albedo, and the asymmetry parameter on the table's wavelength grid
    void interpolateTable(
      const double particle_radius,
      const double sigma_g,
      std::vector<double>& log_cross_section,
      std::vector<double>& single_scattering,
      std::vector<double>& asym_param) const;

    const std::vector<double>& tableRadii() const {return table_radius;}
    const std::vector<double>& tableSigma() const {return table_sigma;}
    const std::vector<double>& tableWavelengths() const {return table_wavelength;}

  protected:
    double reference_wavelength = 0;
    bool fixed_sigma = false;
    double fixed_sigma_value = 1.0;

    //tabulated data; the 3D arrays are stored flat with the wavelength
    //as the fastest index: [(i_r * nb_sigma + i_sigma) * nb_wavelengths + i_wl]
    size_t nb_radii = 0;
    size_t nb_sigma = 0;
    size_t nb_wavelengths = 0;

    std::vector<double> table_radius;
    std::vector<double> table_sigma;
    std::vector<double> table_wavelength;

    std::vector<double> table_log_radius;
    std::vector<double> table_log_wavelength;

    std::vector<double> table_log_cross_section;
    std::vector<double> table_single_scattering;
    std::vector<double> table_asym_param;

    //interpolation position of the reference wavelength in the table
    size_t reference_index = 0;
    double reference_weight = 0;

    //interpolation positions of the spectral grid points in the table;
    //computed on the first call and cached as long as the grid stays the same
    SpectralGrid* mapped_spectral_grid = nullptr;
    size_t mapped_nb_spectral_points = 0;
    std::vector<int> grid_index;
    std::vector<double> grid_weight;

    //GPU copies of the grid mapping and of the per-call interpolated table columns
    int* grid_index_dev = nullptr;
    double* grid_weight_dev = nullptr;
    double* log_cross_section_dev = nullptr;
    double* single_scattering_dev_table = nullptr;
    double* asym_param_dev_table = nullptr;

    void readTable(const std::string& file_path);

    void readParameters(
      const std::vector<double>& parameters,
      double& cloud_optical_depth,
      double& particle_radius,
      double& sigma_g,
      double& cloud_top_pressure,
      double& cloud_bottom_pressure) const;

    void mapSpectralGrid(SpectralGrid* spectral_grid);
    void mapSpectralGridGPU(SpectralGrid* spectral_grid);

    //position of a value in an ascending axis, clamped to the axis range
    static void interpolationPosition(
      const std::vector<double>& axis,
      const double value,
      size_t& index,
      double& weight);

    //log of the extinction cross section interpolated to the reference wavelength
    double referenceLogCrossSection(
      const std::vector<double>& log_cross_section) const;
};


}
#endif 
