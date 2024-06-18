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


#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <omp.h>

#include "discrete_ordinate.h"

#include "../forward_model/atmosphere/atmosphere.h"
#include "../additional/aux_functions.h"
#include "../additional/physical_const.h"
#include "../additional/quadrature.h"
#include "../spectral_grid/spectral_grid.h"
#include "../additional/exceptions.h"

extern "C" {
  #include "cdisort_src/cdisort_macros.h"
  #include "cdisort_src/cdisort.h"
}



namespace helios{

/*
 * Disort-specific shift macros.
 * Using unit-offset shift macros to match Fortran version
 */
#undef  DTAUC
#define DTAUC(lc)  ds.dtauc[lc-1]
#undef  PHI
#define PHI(j)     ds.phi[j-1]
#undef  PMOM
#define PMOM(k,lc) ds.pmom[k+(lc-1)*(ds.nmom_nstr+1)]
#undef  SSALB
#define SSALB(lc)  ds.ssalb[lc-1]
#undef  TEMPER
#define TEMPER(lc) ds.temper[lc]
#undef  UMU
#define UMU(iu)    ds.umu[iu-1]
#undef  UTAU
#define UTAU(lu)   ds.utau[lu-1]


DiscreteOrdinates::DiscreteOrdinates(
  SpectralGrid* spectral_grid_ptr,
  const size_t nb_streams,
  const size_t nb_grid_points,
  const bool use_gpu)
{ 
  if (use_gpu)
  {
    std::string error_message = "Radiative transfer model CDISORT cannot run on the GPU\n";
    throw InvalidInput(std::string ("DiscreteOrdinates::DiscreteOrdinates"), error_message);
  }

  spectral_grid = spectral_grid_ptr;
  initDISORT(nb_streams, nb_grid_points-1);
}



void DiscreteOrdinates::calcSpectrum(
  const Atmosphere& atmosphere,
  const std::vector< std::vector<double> >& absorption_coeff, 
  const std::vector< std::vector<double> >& scattering_coeff,
  const std::vector< std::vector<double> >& cloud_optical_depth,
  const std::vector< std::vector<double> >& cloud_single_scattering,
  const std::vector< std::vector<double> >& cloud_asym_param,
  const double spectrum_scaling,
  std::vector<double>& spectrum)
{ 
  for (size_t i=0; i<ds.size(); ++i)
  {
    setTemperatureStructure(atmosphere.temperature, atmosphere.temperature[0]);
  }


  #pragma omp parallel for
  for (size_t i=0; i<spectrum.size(); ++i)
    spectrum[i] = calcSpectrum(absorption_coeff[i], scattering_coeff[i], cloud_optical_depth[i], atmosphere.altitude, i) * spectrum_scaling;
}




double DiscreteOrdinates::calcSpectrum(
  const std::vector<double> absorption_coeff,
  const std::vector<double> scattering_coeff,
  const std::vector<double>& cloud_optical_depth,
  const std::vector<double>& vertical_grid,
  const size_t nu_index)
{
  size_t nb_grid_points = absorption_coeff.size();

  std::vector<double> optical_depth(nb_grid_points-1, 0.0);


  for (size_t j=0; j<nb_grid_points - 1; ++j)
    optical_depth[j] = (vertical_grid[j+1] - vertical_grid[j]) * (absorption_coeff[j+1] + absorption_coeff[j])/2.;

  if (cloud_optical_depth.size() != 0)
    for (size_t i=0; i<nb_grid_points-1; ++i) optical_depth[i] += cloud_optical_depth[i];


  const double wavenumber = spectral_grid->wavenumber_list[nu_index];
  std::vector<double> single_scattering_albedo(nb_grid_points-1, 0.0);
  double surface_albedo = 0;


  setOpticalDepth(wavenumber, optical_depth, single_scattering_albedo, single_scattering_albedo, surface_albedo);

  double incident_radiation = 0;
  double zenith_angle = 0.5;


  std::vector<double> flux_down(nb_grid_points, 0.0);
  std::vector<double> flux_up(nb_grid_points, 0.0);
  std::vector<double> mean_intensity(nb_grid_points, 0.0);

  calcRadiativeTransfer(incident_radiation, zenith_angle, flux_up, flux_down, mean_intensity);

  return flux_up.back();
}



void DiscreteOrdinates::calcRadiativeTransfer(
  double incident_stellar_radiation,
  double zenith_angle,
  std::vector<double>& flux_up,
  std::vector<double>& flux_down,
  std::vector<double>& mean_intensity)
{ 
  const int thread(omp_get_thread_num());

  ds[thread].bc.fbeam = incident_stellar_radiation;
  ds[thread].bc.umu0  = zenith_angle;


  runDISORT(flux_up, flux_down, mean_intensity);
}



void DiscreteOrdinates::setTemperatureStructure(
  const std::vector<double>& temperature_structure,
  const double& surface_temperature)
{

  for (size_t j=0; j<ds.size(); ++j)
  {
    ds[j].bc.btemp = surface_temperature;

    for (size_t i=0; i<temperature_structure.size(); i++)
      ds[j].temper[i] = temperature_structure[temperature_structure.size() - i - 1];
  }

}



void DiscreteOrdinates::setOpticalDepth(
  const double wavenumber_input,
  const std::vector<double>& optical_depth,
	const std::vector<double>& single_scattering_albedo,
  const std::vector<double>& asymmetry_parameter,
	const double surface_albedo)
{ 
  const int thread(omp_get_thread_num());

  ds[thread].bc.albedo  = surface_albedo;

  ds[thread].wvnmlo = wavenumber_input;
  ds[thread].wvnmhi = wavenumber_input;


  for (int lc = 0; lc < ds[thread].nlyr; lc++)
  {
    ds[thread].dtauc[lc] = optical_depth[ds[thread].nlyr - lc - 1];
    ds[thread].ssalb[lc] = single_scattering_albedo[ds[thread].nlyr - lc - 1];
  }

  for (int lc = 0; lc < ds[thread].nlyr; lc++)
  {
    double gg = asymmetry_parameter[ds[thread].nlyr - lc - 1];

    if (gg > 0)
      c_getmom(HENYEY_GREENSTEIN,gg,ds[thread].nmom,&ds[thread].pmom[0 + (lc)*(ds[thread].nmom_nstr+1)]);
    else
      c_getmom(RAYLEIGH,gg,ds[thread].nmom,&ds[thread].pmom[0 + (lc)*(ds[thread].nmom_nstr+1)]);
  }

}


void DiscreteOrdinates::initDISORT(unsigned int nb_streams, unsigned int nb_layers)
{ 
  disort_state ds_gen;

  ds_gen.nstr   = nb_streams;
  ds_gen.nphase = ds_gen.nstr;
  ds_gen.nlyr   = nb_layers;
  ds_gen.nmom   = nb_streams;
  ds_gen.ntau   = 0;
  ds_gen.numu   = 0;
  ds_gen.nphi   = 0;
  ds_gen.flag.usrtau = 0;
  ds_gen.flag.planck = 1;


  ds_gen.accur = 0.;
  //ds_gen.flag.prnt[0]=TRUE, ds_gen.flag.prnt[1]=TRUE, ds_gen.flag.prnt[2]=TRUE, ds_gen.flag.prnt[3]=TRUE, ds_gen.flag.prnt[4]=TRUE;
  ds_gen.flag.prnt[0]=FALSE, ds_gen.flag.prnt[1]=FALSE, ds_gen.flag.prnt[2]=FALSE, ds_gen.flag.prnt[3]=FALSE, ds_gen.flag.prnt[4]=FALSE;

  ds_gen.flag.ibcnd  = GENERAL_BC;
  ds_gen.flag.usrtau = FALSE;
  ds_gen.flag.usrang = FALSE;
  ds_gen.flag.lamber = TRUE;
  ds_gen.flag.onlyfl = TRUE;
  ds_gen.flag.quiet  = TRUE;
  ds_gen.flag.spher  = FALSE;
  ds_gen.flag.general_source = FALSE;
  ds_gen.flag.output_uum = FALSE;
  ds_gen.flag.intensity_correction = TRUE;
  ds_gen.flag.old_intensity_correction = FALSE;


  ds_gen.bc.fisot = 0;
  ds_gen.bc.phi0  = 0.0;
  ds_gen.bc.fluor = 0.;

  ds_gen.flag.brdf_type = BRDF_NONE;

  ds_gen.bc.ttemp   = 0;
  ds_gen.bc.temis   = 0;


  const int nb_cores(omp_get_max_threads());
  
  ds.assign(nb_cores, ds_gen);
  out.resize(nb_cores);


  /* Allocate memory */
  for (int i=0; i<nb_cores; ++i)
  {
    c_disort_state_alloc(&ds[i]);
    c_disort_out_alloc(&ds[i],&out[i]);
  }

}



void DiscreteOrdinates::finaliseDISORT()
{
   const int nb_cores(omp_get_max_threads());

  //free allocated memory
  for (int i=0; i<nb_cores; ++i)
  {
    c_disort_out_free(&ds[i],&out[i]);
    c_disort_state_free(&ds[i]);
  }

}



void DiscreteOrdinates::runDISORT(
  std::vector<double>& flux_up,
  std::vector<double>& flux_down,
  std::vector<double>& mean_intensity)
{ 
  const int thread(omp_get_thread_num());

  if (c_planck_func2(ds[thread].wvnmlo, ds[thread].wvnmlo, ds[thread].bc.btemp) > 1.e-35)
    ds[thread].flag.planck = TRUE;
  else
    ds[thread].flag.planck = FALSE;


  c_disort(&ds[thread],&out[thread]);


  for (int lc = 0; lc <= ds[thread].nlyr; lc++)
  {
    flux_down[ds[thread].nlyr - lc] = out[thread].rad[lc].rfldir + out[thread].rad[lc].rfldn;
    flux_up[ds[thread].nlyr - lc] = out[thread].rad[lc].flup;
    mean_intensity[ds[thread].nlyr - lc] = out[thread].rad[lc].uavg;
  }

}


#undef PMOM
}
