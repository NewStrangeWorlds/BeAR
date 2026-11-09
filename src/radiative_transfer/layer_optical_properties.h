// SPDX-FileCopyrightText: 2024-2026 Daniel Kitzmann
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BeAR (https://github.com/newstrangeworlds/BeAR).

#ifndef _layer_optical_properties_h
#define _layer_optical_properties_h

#include <vector>
#include <cstddef>


namespace bear {


//combined optical properties of a single atmospheric layer
struct LayerOpticalProperties {
  double optical_depth = 0.0;             //gas absorption + Rayleigh scattering + cloud extinction
  double single_scattering = 0.0;         //total single scattering albedo (gas + cloud)
  double cloud_scattering_fraction = 0.0; //fraction of the scattering optical depth due to clouds
  double cloud_asym_param = 0.0;          //asymmetry parameter of the cloud scattering
};



//assembles the optical properties of layer j from the level-based gas coefficients
//and the layer-based cloud properties of a single spectral point
//gas scattering (Rayleigh) is assumed to be isotropic in the sense of the asymmetry
//parameter, the effective asymmetry parameter of the layer therefore is
//cloud_scattering_fraction * cloud_asym_param
inline LayerOpticalProperties layerOpticalProperties(
  const std::vector<double>& altitude,
  const std::vector<double>& absorption_coeff,
  const std::vector<double>& scattering_coeff,
  const std::vector<double>& cloud_optical_depth,
  const std::vector<double>& cloud_single_scattering,
  const std::vector<double>& cloud_asym_param,
  const size_t j)
{
  LayerOpticalProperties p;

  const double dz = altitude[j+1] - altitude[j];

  const double absorption_depth = dz * (absorption_coeff[j+1] + absorption_coeff[j]) / 2.;
  const double gas_scattering_depth = dz * (scattering_coeff[j+1] + scattering_coeff[j]) / 2.;

  double cloud_scattering_depth = 0.0;

  p.optical_depth = absorption_depth + gas_scattering_depth;

  if (cloud_optical_depth.size() != 0)
  {
    p.optical_depth += cloud_optical_depth[j];

    if (cloud_single_scattering.size() != 0)
      cloud_scattering_depth = cloud_optical_depth[j] * cloud_single_scattering[j];

    if (cloud_asym_param.size() != 0)
      p.cloud_asym_param = cloud_asym_param[j];
  }

  if (p.optical_depth < 0.0) p.optical_depth = 0.0;

  const double scattering_depth = gas_scattering_depth + cloud_scattering_depth;

  if (p.optical_depth > 0.0)
    p.single_scattering = scattering_depth / p.optical_depth;

  if (p.single_scattering > 1.0) p.single_scattering = 1.0;

  if (scattering_depth > 0.0)
    p.cloud_scattering_fraction = cloud_scattering_depth / scattering_depth;

  return p;
}


}

#endif
