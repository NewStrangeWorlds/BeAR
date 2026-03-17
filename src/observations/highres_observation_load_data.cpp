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


#include "highres_observation.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

#include "../additional/exceptions.h"


namespace bear {


void HighResObservation::loadDataFile(const std::string& file_path)
{
  std::fstream file(file_path.c_str(), std::ios::in);

  if (file.fail())
  {
    std::string error_message =
      "Couldn't open high-res observation file: " + file_path + "\n";
    throw InvalidInput(std::string("HighResObservation::loadDataFile"), error_message);
  }

  std::string line;

  // Read header fields
  while (std::getline(file, line))
  {
    if (line.empty() || line[0] != '#') continue;

    if (line == "#name")
    {
      std::getline(file, line);
      observation_name = line;
    }
    else if (line == "#resolving_power")
    {
      file >> resolving_power;
      std::getline(file, line);  // consume newline
    }
    else if (line == "#nb_orders")
    {
      file >> nb_orders;
      std::getline(file, line);
    }
    else if (line == "#nb_exposures")
    {
      file >> nb_exposures;
      std::getline(file, line);
    }
    else if (line == "#orbital_phases")
    {
      orbital_phases.resize(nb_exposures);
      for (size_t i = 0; i < nb_exposures; ++i)
        file >> orbital_phases[i];
      std::getline(file, line);
    }
    else if (line.substr(0, 6) == "#order")
    {
      // Parse order index (not strictly needed, read sequentially)
      SpectralOrder order;

      // Read wavelengths
      std::getline(file, line);  // should be "#wavelengths"

      std::getline(file, line);
      std::istringstream wl_stream(line);
      double val;
      while (wl_stream >> val)
        order.wavelengths.push_back(val);

      order.nb_pixels = order.wavelengths.size();

      // Read flux
      std::getline(file, line);  // should be "#flux"

      order.flux.resize(nb_exposures);
      for (size_t e = 0; e < nb_exposures; ++e)
      {
        std::getline(file, line);
        std::istringstream flux_stream(line);
        order.flux[e].resize(order.nb_pixels);
        for (size_t p = 0; p < order.nb_pixels; ++p)
          flux_stream >> order.flux[e][p];
      }

      // Check for optional #flux_uncertainties section
      std::streampos pos = file.tellg();
      std::getline(file, line);

      if (line == "#flux_uncertainties")
      {
        has_flux_uncertainties = true;

        order.flux_uncertainties.resize(nb_exposures);
        for (size_t e = 0; e < nb_exposures; ++e)
        {
          std::getline(file, line);
          std::istringstream unc_stream(line);
          order.flux_uncertainties[e].resize(order.nb_pixels);
          for (size_t p = 0; p < order.nb_pixels; ++p)
          {
            unc_stream >> order.flux_uncertainties[e][p];
            // Floor to prevent division by zero
            if (order.flux_uncertainties[e][p] < 1e-30)
              order.flux_uncertainties[e][p] = 1e-30;
          }
        }
      }
      else
      {
        // No uncertainties — seek back
        file.seekg(pos);
      }

      spectral_orders.push_back(std::move(order));
    }
    else if (line == "#filtering_basis")
    {
      has_filtering = true;
      loadFilteringBasis(file_path);
      break;  // filtering section must be last
    }
  }

  file.close();

  if (spectral_orders.size() != nb_orders)
  {
    std::string error_message =
      "Expected " + std::to_string(nb_orders)
      + " orders but found " + std::to_string(spectral_orders.size())
      + " in file: " + file_path + "\n";
    throw InvalidInput(std::string("HighResObservation::loadDataFile"), error_message);
  }
}


void HighResObservation::loadFilteringBasis(const std::string& file_path)
{
  std::fstream file(file_path.c_str(), std::ios::in);

  if (file.fail())
    throw InvalidInput(
      std::string("HighResObservation::loadFilteringBasis"),
      "Couldn't open file: " + file_path + "\n");

  std::string line;

  // Skip to #filtering_basis section
  while (std::getline(file, line))
  {
    if (line == "#filtering_basis")
      break;
  }

  // Read number of basis vectors
  while (std::getline(file, line))
  {
    if (line == "#nb_basis_vectors")
    {
      file >> nb_basis_vectors;
      std::getline(file, line);  // consume newline
      break;
    }
  }

  basis_vectors_raw.resize(nb_orders);
  uncertainties.resize(nb_orders);

  size_t orders_read = 0;

  while (std::getline(file, line))
  {
    if (line.empty() || line[0] != '#') continue;

    if (line.substr(0, 6) == "#order")
    {
      if (orders_read >= nb_orders)
        throw InvalidInput(
          std::string("HighResObservation::loadFilteringBasis"),
          "Too many orders in filtering_basis section\n");

      size_t ord = orders_read;

      // Read basis vectors: nb_exposures lines, each with nb_basis_vectors values
      std::getline(file, line);  // should be "#basis_vectors"

      basis_vectors_raw[ord].resize(nb_exposures * nb_basis_vectors);

      for (size_t e = 0; e < nb_exposures; ++e)
      {
        std::getline(file, line);
        std::istringstream stream(line);

        for (size_t b = 0; b < nb_basis_vectors; ++b)
          stream >> basis_vectors_raw[ord][e * nb_basis_vectors + b];
      }

      // Check for optional uncertainties
      std::streampos pos = file.tellg();
      std::getline(file, line);

      if (line == "#uncertainties")
      {
        uncertainties[ord].resize(nb_exposures);
        std::getline(file, line);
        std::istringstream unc_stream(line);

        for (size_t e = 0; e < nb_exposures; ++e)
          unc_stream >> uncertainties[ord][e];
      }
      else
      {
        // No uncertainties — seek back
        file.seekg(pos);
        uncertainties[ord].clear();
      }

      ++orders_read;
    }
  }

  file.close();

  if (orders_read != nb_orders)
    throw InvalidInput(
      std::string("HighResObservation::loadFilteringBasis"),
      "Expected " + std::to_string(nb_orders) + " orders in filtering_basis "
      "but found " + std::to_string(orders_read) + "\n");

  std::cout << "  Filtering basis loaded: " << nb_basis_vectors
            << " basis vectors per order\n";
}


}
