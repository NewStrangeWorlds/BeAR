
#include "species_definition.h"
#include "opacity_species.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <assert.h>
#include <omp.h>


#include "../config/global_config.h"
#include "../spectral_grid/spectral_grid.h"


namespace helios{


bool GasH2::calcRalyleighCrossSections(std::vector<double>& cross_sections)
{
  unsigned int nb_spectral_points = spectral_grid->nbSpectralPoints();


  #pragma omp parallel for
  for (unsigned int j=0; j<nb_spectral_points; j++)
  {
    double king_correction_factor = 1.0;

    double refractive_index = (13.58e-5 * (1. + 7.52e-3/std::pow(spectral_grid->wavelength_list[j],2))) + 1.;


    double reference_density = 2.651629e19; //molecules cm^-3

    cross_sections[j] = generalRayleighCrossSection(reference_density, refractive_index, king_correction_factor, spectral_grid->wavenumber_list[j]);
  }

  return true;
}



bool GasHe::calcRalyleighCrossSections(std::vector<double>& cross_sections)
{
  unsigned int nb_spectral_points = spectral_grid->nbSpectralPoints();


  #pragma omp parallel for
  for (unsigned int j=0; j<nb_spectral_points; j++)
  {
    double king_correction_factor = 1.0;

    double refractive_index = (2283. + 1.8102e13 / (1.5342e10 - spectral_grid->wavenumber_list[j]*spectral_grid->wavenumber_list[j])) * 1e-8 + 1;


    double reference_density = 2.546899e19; //molecules cm^-3

    cross_sections[j] = generalRayleighCrossSection(reference_density, refractive_index, king_correction_factor, spectral_grid->wavenumber_list[j]);
  }

  return true;
}


bool GasH::calcRalyleighCrossSections(std::vector<double>& cross_sections)
{
  unsigned int nb_spectral_points = spectral_grid->nbSpectralPoints();


  #pragma omp parallel for
  for (unsigned int i=0; i<nb_spectral_points; i++)
  {
    const double sigma_thomson = 0.665e-24; //Thomson scattering cross-section in cm2
    const double lambda_lyman = 0.0912; //wavelength of the Lyman limit in micron

    const double lambda_fraction = lambda_lyman / spectral_grid->wavelength_list[i];
    
    cross_sections[i] = 8.41e-25 * std::pow(lambda_fraction, 4) + 3.37e-24 * std::pow(lambda_fraction, 6) + 4.71e-22 * std::pow(lambda_fraction, 14); //in cm2
  }

  return true;
}


}