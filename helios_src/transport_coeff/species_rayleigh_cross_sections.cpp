
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


bool GasCO::calcRalyleighCrossSections(std::vector<double>& cross_sections)
{
  unsigned int nb_spectral_points = spectral_grid->nbSpectralPoints();
  
  double reference_density = 2.546899e19; //molecules cm^-3
  double king_correction_factor = 1.0;

  #pragma omp parallel for
  for (unsigned int j=0; j<nb_spectral_points; j++)
  {
    if (spectral_grid->wavelength_list[j] > 2.0) continue;

    double refractive_index = (22851. + 0.456e14 / (71427.0*71427.0 - spectral_grid->wavenumber_list[j]*spectral_grid->wavenumber_list[j])) * 1e-8 + 1;

    cross_sections[j] = generalRayleighCrossSection(reference_density, refractive_index, king_correction_factor, spectral_grid->wavenumber_list[j]);
  }


  return true;
}



bool GasCO2::calcRalyleighCrossSections(std::vector<double>& cross_sections)
{
  unsigned int nb_spectral_points = spectral_grid->nbSpectralPoints();
  
  double reference_density = 2.546899e19; //molecules cm^-3

  #pragma omp parallel for
  for (unsigned int j=0; j<nb_spectral_points; j++)
  {
    if (spectral_grid->wavelength_list[j] > 2.0) continue;

    double king_correction_factor = 1.1364 + 25.3e-12 * spectral_grid->wavenumber_list[j]*spectral_grid->wavenumber_list[j];

    double refractive_index = (5799.25 / (128908.9*128908.9 - spectral_grid->wavenumber_list[j]*spectral_grid->wavenumber_list[j]) 
                              + 120.05 / (89223.8*89223.8 - spectral_grid->wavenumber_list[j]*spectral_grid->wavenumber_list[j]) 
                              + 5.3334 / (75037.5*75037.5 - spectral_grid->wavenumber_list[j]*spectral_grid->wavenumber_list[j]) 
                              + 4.3244 / (67837.7*67837.7 - spectral_grid->wavenumber_list[j]*spectral_grid->wavenumber_list[j]) 
                              + 0.1218145e-6 / (2418.136*2418.136 - spectral_grid->wavenumber_list[j]*spectral_grid->wavenumber_list[j]))
                              * 1.1427e3 + 1.0;

    cross_sections[j] = generalRayleighCrossSection(reference_density, refractive_index, king_correction_factor, spectral_grid->wavenumber_list[j]);
  }

  return true;
}



bool GasCH4::calcRalyleighCrossSections(std::vector<double>& cross_sections)
{
  unsigned int nb_spectral_points = spectral_grid->nbSpectralPoints();
  
  double reference_density = 2.546899e19; //molecules cm^-3
  double king_correction_factor = 1.0;

  #pragma omp parallel for
  for (unsigned int j=0; j<nb_spectral_points; j++)
  {
    if (spectral_grid->wavelength_list[j] > 2.0) continue;

    double refractive_index = (46662. + 4.02e-6 * spectral_grid->wavenumber_list[j]*spectral_grid->wavenumber_list[j]) * 1e-8 + 1;

    cross_sections[j] = generalRayleighCrossSection(reference_density, refractive_index, king_correction_factor, spectral_grid->wavenumber_list[j]);
  }


  return true;
}



bool GasH2O::calcRalyleighCrossSections(std::vector<double>& cross_sections)
{
  unsigned int nb_spectral_points = spectral_grid->nbSpectralPoints();

  //this is the number density of water at standard temperature and pressure
  //determined from the Avogradro constant and the properties of water at STP
  double reference_density = 3.34279671749673e+22;
  
  //values for water at STP
  double delta = 1.0;
  double theta = 1.0;

  double lambda_uv = 0.229202;
  double lambda_ir = 5.432937;

  std::vector<double> a_coeff = {0.244257733, 0.974634476e-2, -0.373234996e-2, 0.268678472e-3, 0.158920570e-2, 0.245934259e-2, 0.900704920, -0.166626219e-1};
  
  double king_correction_factor = (6 + 3 * 3e-4) / (6 - 7 * 3e-4);


  #pragma omp parallel for
  for (unsigned int j=0; j<nb_spectral_points; j++)
  {
    if (spectral_grid->wavelength_list[j] > 2.0) continue;

    double lambda = spectral_grid->wavelength_list[j] / 0.589;

    double a_factor = delta * (a_coeff[0] 
                             + a_coeff[1]*delta 
                             + a_coeff[2]*theta 
                             + a_coeff[3]*lambda*lambda*theta 
                             + a_coeff[4]*std::pow(lambda,-2) 
                             + a_coeff[5] / (lambda*lambda - lambda_uv*lambda_uv) 
                             + a_coeff[6] / (lambda*lambda - lambda_ir*lambda_ir) + a_coeff[7]*delta*delta);
  
    double refractive_index = std::pow(((2 * a_factor + 1)/(1 - a_factor)), 0.5);


    cross_sections[j] = generalRayleighCrossSection(reference_density, refractive_index, king_correction_factor, spectral_grid->wavenumber_list[j]);
  }


  return true;
}



}