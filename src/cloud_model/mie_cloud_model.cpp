// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#include "mie_cloud_model.h"

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdlib>

#include "grey_cloud_model.h"
#include "../forward_model/atmosphere/atmosphere.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../additional/exceptions.h"


namespace bear{


//model options in the config file:
//  parameters[0]   table file (relative to the retrieval folder)
//  parameters[1]   reference wavelength in micron for the optical depth parameter
//  "fb"            fixed cloud bottom, one scale height below the cloud top
//  "sigma_g=<x>"   fixes the width of the size distribution instead of retrieving it
MieCloudModel::MieCloudModel(
  const std::vector<std::string>& parameters,
  const std::string& retrieval_folder)
 : GreyCloudModel(std::vector<std::string>{})
{
  if (parameters.size() < 2)
  {
    std::string error_message =
      "The Mie cloud model requires at least two parameters (table file, reference wavelength), but found "
      + std::to_string(parameters.size()) + "\n";
    throw InvalidInput(std::string ("forward_model.toml"), error_message);
  }

  const std::string table_file = retrieval_folder + parameters[0];

  reference_wavelength = std::stod(parameters[1]);

  if (reference_wavelength <= 0)
  {
    std::string error_message =
      "The reference wavelength of the Mie cloud model must be positive, found " + parameters[1] + "\n";
    throw InvalidInput(std::string ("forward_model.toml"), error_message);
  }

  for (size_t i=2; i<parameters.size(); ++i)
  {
    const std::string& p = parameters[i];

    if (p == "fb")
      fixed_bottom = true;
    else if (p.rfind("sigma_g=", 0) == 0)
    {
      fixed_sigma = true;
      fixed_sigma_value = std::stod(p.substr(8));
    }
    else
    {
      std::string error_message =
        "Unknown option \"" + p + "\" for the Mie cloud model (allowed: fb, sigma_g=<value>)\n";
      throw InvalidInput(std::string ("forward_model.toml"), error_message);
    }
  }


  readTable(table_file);


  if (reference_wavelength < table_wavelength.front() || reference_wavelength > table_wavelength.back())
  {
    std::string error_message =
      "The reference wavelength " + parameters[1] + " of the Mie cloud model is outside the wavelength range of the table "
      + table_file + " (" + std::to_string(table_wavelength.front()) + " - " 
      + std::to_string(table_wavelength.back()) + " micron)\n";
    throw InvalidInput(std::string ("forward_model.toml"), error_message);
  }

  interpolationPosition(table_log_wavelength, std::log(reference_wavelength), reference_index, reference_weight);


  if (fixed_sigma && (fixed_sigma_value < table_sigma.front() || fixed_sigma_value > table_sigma.back()))
  {
    std::string error_message =
      "The fixed sigma_g value " + std::to_string(fixed_sigma_value) + " of the Mie cloud model is outside the range of the table "
      + table_file + " (" + std::to_string(table_sigma.front()) + " - " 
      + std::to_string(table_sigma.back()) + ")\n";
    throw InvalidInput(std::string ("forward_model.toml"), error_message);
  }


  //populate the parameter names in the same order readParameters() reads them
  parameter_names = {
    "cloud_optical_depth",
    "cloud_particle_size"};

  if (fixed_sigma == false)
    parameter_names.push_back("cloud_sigma_g");

  parameter_names.push_back("cloud_top_pressure");

  if (fixed_bottom == false)
    parameter_names.push_back("cloud_bottom_fraction");


  std::cout << "Mie cloud model: read table " << table_file << " with " 
            << nb_radii << " particle radii (" << table_radius.front() << " - " << table_radius.back() << " micron), "
            << nb_sigma << " sigma_g values (" << table_sigma.front() << " - " << table_sigma.back() << "), and "
            << nb_wavelengths << " wavelengths (" << table_wavelength.front() << " - " << table_wavelength.back() << " micron)\n";
}



MieCloudModel::~MieCloudModel()
{
  if (grid_index_dev != nullptr) deleteFromDevice(grid_index_dev);
  if (grid_weight_dev != nullptr) deleteFromDevice(grid_weight_dev);
  if (log_cross_section_dev != nullptr) deleteFromDevice(log_cross_section_dev);
  if (single_scattering_dev_table != nullptr) deleteFromDevice(single_scattering_dev_table);
  if (asym_param_dev_table != nullptr) deleteFromDevice(asym_param_dev_table);
}



//reads the ASCII table written by tools/mie_table_to_ascii.py
//lines starting with '#' are treated as comments
void MieCloudModel::readTable(const std::string& file_path)
{
  std::fstream file(file_path.c_str(), std::ios::in);

  if (file.fail())
    throw FileNotFound(std::string ("MieCloudModel::readTable"), file_path);

  //reads the next non-comment token from the file
  auto nextToken = [&file](std::string& token) -> bool {
    while (file >> token)
    {
      if (token[0] == '#')
      {
        std::string rest_of_line;
        std::getline(file, rest_of_line);
        continue;
      }

      return true;
    }

    return false;
  };

  auto nextValue = [&nextToken, &file_path](double& value) {
    std::string token;

    if (!nextToken(token))
      throw InvalidInput(std::string ("MieCloudModel::readTable"), 
        "Unexpected end of the Mie cloud table file " + file_path + "\n");

    char* end = nullptr;
    value = std::strtod(token.c_str(), &end);

    if (end == token.c_str() || *end != '\0')
      throw InvalidInput(std::string ("MieCloudModel::readTable"), 
        "Could not parse the value \"" + token + "\" in the Mie cloud table file " + file_path + "\n");
  };

  auto readAxis = [&nextValue, &file_path](std::vector<double>& axis, const std::string& name) {
    for (auto & a : axis) nextValue(a);

    for (size_t i=1; i<axis.size(); ++i)
      if (axis[i] <= axis[i-1])
        throw InvalidInput(std::string ("MieCloudModel::readTable"), 
          "The " + name + " axis in the Mie cloud table file " + file_path + " is not strictly ascending\n");
  };


  double value = 0;

  nextValue(value); nb_radii = static_cast<size_t>(value);
  nextValue(value); nb_sigma = static_cast<size_t>(value);
  nextValue(value); nb_wavelengths = static_cast<size_t>(value);

  if (nb_radii < 1 || nb_sigma < 1 || nb_wavelengths < 2)
    throw InvalidInput(std::string ("MieCloudModel::readTable"), 
      "Invalid table dimensions in the Mie cloud table file " + file_path + "\n");


  table_radius.resize(nb_radii);
  table_sigma.resize(nb_sigma);
  table_wavelength.resize(nb_wavelengths);

  readAxis(table_radius, "particle radius");
  readAxis(table_sigma, "sigma_g");
  readAxis(table_wavelength, "wavelength");

  if (table_radius.front() <= 0)
    throw InvalidInput(std::string ("MieCloudModel::readTable"), 
      "The particle radii in the Mie cloud table file " + file_path + " must be positive\n");


  const size_t nb_entries = nb_radii * nb_sigma * nb_wavelengths;

  table_log_cross_section.resize(nb_entries);
  table_single_scattering.resize(nb_entries);
  table_asym_param.resize(nb_entries);

  for (size_t i=0; i<nb_radii; ++i)
    for (size_t j=0; j<nb_sigma; ++j)
    {
      const size_t offset = (i*nb_sigma + j) * nb_wavelengths;

      for (size_t k=0; k<nb_wavelengths; ++k)
      {
        nextValue(value);

        if (value <= 0)
          throw InvalidInput(std::string ("MieCloudModel::readTable"), 
            "Non-positive extinction cross section in the Mie cloud table file " + file_path + "\n");

        table_log_cross_section[offset + k] = std::log(value);
      }

      for (size_t k=0; k<nb_wavelengths; ++k)
      {
        nextValue(value);
        table_single_scattering[offset + k] = std::min(std::max(value, 0.0), 1.0);
      }

      for (size_t k=0; k<nb_wavelengths; ++k)
      {
        nextValue(value);
        table_asym_param[offset + k] = std::min(std::max(value, -1.0), 1.0);
      }
    }

  file.close();


  table_log_radius.resize(nb_radii);
  table_log_wavelength.resize(nb_wavelengths);

  for (size_t i=0; i<nb_radii; ++i)
    table_log_radius[i] = std::log(table_radius[i]);

  for (size_t i=0; i<nb_wavelengths; ++i)
    table_log_wavelength[i] = std::log(table_wavelength[i]);
}



//finds the interval of an ascending axis that contains a value and the
//linear interpolation weight of the upper interval point
//values outside the axis are clamped to the axis boundaries
void MieCloudModel::interpolationPosition(
  const std::vector<double>& axis,
  const double value,
  size_t& index,
  double& weight)
{
  if (axis.size() < 2 || value <= axis.front())
  {
    index = 0;
    weight = 0;
    return;
  }

  if (value >= axis.back())
  {
    index = axis.size() - 2;
    weight = 1.0;
    return;
  }

  auto it = std::upper_bound(axis.begin(), axis.end(), value);

  index = std::distance(axis.begin(), it) - 1;
  weight = (value - axis[index]) / (axis[index+1] - axis[index]);
}



//bilinear interpolation in log(r_g) and sigma_g
//the extinction cross section is interpolated in log space
void MieCloudModel::interpolateTable(
  const double particle_radius,
  const double sigma_g,
  std::vector<double>& log_cross_section,
  std::vector<double>& single_scattering,
  std::vector<double>& asym_param) const
{
  size_t ir = 0;
  double wr = 0;
  interpolationPosition(table_log_radius, std::log(particle_radius), ir, wr);

  size_t is = 0;
  double ws = 0;
  interpolationPosition(table_sigma, sigma_g, is, ws);

  //a single-entry axis only has one interpolation point
  const size_t ir1 = (nb_radii > 1) ? ir + 1 : ir;
  const size_t is1 = (nb_sigma > 1) ? is + 1 : is;

  const double w00 = (1.0 - wr) * (1.0 - ws);
  const double w01 = (1.0 - wr) * ws;
  const double w10 = wr * (1.0 - ws);
  const double w11 = wr * ws;

  const size_t o00 = (ir  * nb_sigma + is ) * nb_wavelengths;
  const size_t o01 = (ir  * nb_sigma + is1) * nb_wavelengths;
  const size_t o10 = (ir1 * nb_sigma + is ) * nb_wavelengths;
  const size_t o11 = (ir1 * nb_sigma + is1) * nb_wavelengths;

  log_cross_section.resize(nb_wavelengths);
  single_scattering.resize(nb_wavelengths);
  asym_param.resize(nb_wavelengths);

  for (size_t k=0; k<nb_wavelengths; ++k)
  {
    log_cross_section[k] = 
        w00 * table_log_cross_section[o00 + k] + w01 * table_log_cross_section[o01 + k]
      + w10 * table_log_cross_section[o10 + k] + w11 * table_log_cross_section[o11 + k];

    single_scattering[k] = 
        w00 * table_single_scattering[o00 + k] + w01 * table_single_scattering[o01 + k]
      + w10 * table_single_scattering[o10 + k] + w11 * table_single_scattering[o11 + k];

    asym_param[k] = 
        w00 * table_asym_param[o00 + k] + w01 * table_asym_param[o01 + k]
      + w10 * table_asym_param[o10 + k] + w11 * table_asym_param[o11 + k];
  }
}



double MieCloudModel::referenceLogCrossSection(
  const std::vector<double>& log_cross_section) const
{
  return (1.0 - reference_weight) * log_cross_section[reference_index] 
         + reference_weight * log_cross_section[reference_index + 1];
}



//computes the interpolation positions of the spectral grid points in the
//table's wavelength grid; points outside the table are clamped to its edges
void MieCloudModel::mapSpectralGrid(SpectralGrid* spectral_grid)
{
  const size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  if (mapped_spectral_grid == spectral_grid && mapped_nb_spectral_points == nb_spectral_points)
    return;

  grid_index.assign(nb_spectral_points, 0);
  grid_weight.assign(nb_spectral_points, 0.0);

  size_t nb_outside = 0;

  for (size_t i=0; i<nb_spectral_points; ++i)
  {
    const double wavelength = spectral_grid->wavelength_list[i];

    if (wavelength < table_wavelength.front() || wavelength > table_wavelength.back())
      nb_outside++;

    size_t index = 0;
    double weight = 0;
    interpolationPosition(table_log_wavelength, std::log(wavelength), index, weight);

    grid_index[i] = static_cast<int>(index);
    grid_weight[i] = weight;
  }

  if (nb_outside > 0)
    std::cout << "Warning: Mie cloud model: " << nb_outside << " of " << nb_spectral_points 
              << " spectral grid points lie outside the wavelength range of the table (" 
              << table_wavelength.front() << " - " << table_wavelength.back() 
              << " micron); the table edge values are used there.\n";

  mapped_spectral_grid = spectral_grid;
  mapped_nb_spectral_points = nb_spectral_points;

  //invalidate a potential GPU copy of a previous mapping
  if (grid_index_dev != nullptr) deleteFromDevice(grid_index_dev);
  if (grid_weight_dev != nullptr) deleteFromDevice(grid_weight_dev);
}



void MieCloudModel::readParameters(
  const std::vector<double>& parameters,
  double& cloud_optical_depth,
  double& particle_radius,
  double& sigma_g,
  double& cloud_top_pressure,
  double& cloud_bottom_pressure) const
{
  size_t idx = 0;

  cloud_optical_depth = parameters[idx++];
  particle_radius = parameters[idx++];

  sigma_g = fixed_sigma_value;

  if (fixed_sigma == false)
    sigma_g = parameters[idx++];

  cloud_top_pressure = parameters[idx++];

  cloud_bottom_pressure = 0;

  if (fixed_bottom == false)
    cloud_bottom_pressure = cloud_top_pressure * parameters[idx++];

  if (cloud_optical_depth < 0) cloud_optical_depth = 0;

  //values outside the table are clamped to the table boundaries by the interpolation
}



//calculates the vertical distribution of the cloud layer
//the optical depth at the reference wavelength is distributed evenly over the
//layers between the cloud's top and bottom; its wavelength dependence follows
//the tabulated extinction cross section, while the single scattering albedo and
//the asymmetry parameter are taken directly from the table
void MieCloudModel::opticalProperties(
  const std::vector<double>& parameters, 
  const Atmosphere& atmosphere,
  SpectralGrid* spectral_grid,
  std::vector<std::vector<double>>& optical_depth, 
  std::vector<std::vector<double>>& single_scattering, 
  std::vector<std::vector<double>>& asym_param)
{
  double cloud_optical_depth = 0;
  double particle_radius = 0;
  double sigma_g = 0;
  double cloud_top_pressure = 0;
  double cloud_bottom_pressure = 0;

  readParameters(
    parameters,
    cloud_optical_depth,
    particle_radius,
    sigma_g,
    cloud_top_pressure,
    cloud_bottom_pressure);


  const size_t nb_spectral_points = spectral_grid->nbSpectralPoints();
  const size_t nb_grid_points = atmosphere.nb_grid_points;

  //pre-allocate the data
  optical_depth.assign(nb_spectral_points, std::vector<double>(nb_grid_points-1, 0.0));
  single_scattering.assign(nb_spectral_points, std::vector<double>(nb_grid_points-1, 0.0));
  asym_param.assign(nb_spectral_points, std::vector<double>(nb_grid_points-1, 0.0));


  unsigned int cloud_top_index = 0;
  unsigned int cloud_bottom_index = 0;

  if (fixed_bottom == true)
    cloudPosition(atmosphere, cloud_top_pressure, cloud_top_index, cloud_bottom_index);
  else
    cloudPosition(atmosphere, cloud_top_pressure, cloud_bottom_pressure, cloud_top_index, cloud_bottom_index);


  const double optical_depth_layer_reference = 
    cloud_optical_depth/static_cast<double>(cloud_top_index - cloud_bottom_index);


  std::vector<double> log_cross_section;
  std::vector<double> table_ssa;
  std::vector<double> table_g;

  interpolateTable(particle_radius, sigma_g, log_cross_section, table_ssa, table_g);

  const double log_cross_section_ref = referenceLogCrossSection(log_cross_section);

  mapSpectralGrid(spectral_grid);


  for (size_t j=0; j<nb_spectral_points; ++j)
  {
    const size_t k = grid_index[j];
    const double w = grid_weight[j];

    const double log_cross_section_wl = (1.0 - w) * log_cross_section[k] + w * log_cross_section[k+1];
    const double ssa_wl = (1.0 - w) * table_ssa[k] + w * table_ssa[k+1];
    const double g_wl = (1.0 - w) * table_g[k] + w * table_g[k+1];

    const double optical_depth_layer = 
      optical_depth_layer_reference * std::exp(log_cross_section_wl - log_cross_section_ref);

    for (size_t i=cloud_top_index; i>cloud_bottom_index; --i)
    {
      optical_depth[j][i] = optical_depth_layer;
      single_scattering[j][i] = ssa_wl;
      asym_param[j][i] = g_wl;
    }
  }
}


}
