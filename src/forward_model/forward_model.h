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


#ifndef _forward_model_h
#define _forward_model_h

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "../additional/exceptions.h"
#include "../retrieval/priors.h"
#include "../config/global_config.h"
#include "../spectral_grid/spectral_grid.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../observations/observations.h"


namespace bear {


//abstract class for the forward model
//a derived class *has to* implement all the various, virtual methods
class ForwardModel{
  public:
  ForwardModel (
      GlobalConfig* config_, 
      SpectralGrid* spectral_grid_,
      std::vector<Observation>& observations_);
    virtual ~ForwardModel() {}
    //calculate a model on the CPU
    //the return value signals the retrieval to neglect this model
    virtual bool calcModel(
      const std::vector<double>& parameter, 
      std::vector<double>& spectrum, 
      std::vector<std::vector<double>>& spectrum_obs) = 0;
    //calculate a model on the GPU
    //the return value signals the retrieval to neglect this model
    virtual bool calcModelGPU(
      const std::vector<double>& parameter, 
      double* spectrum, 
      std::vector<double*>& spectrum_obs) = 0;
    //model-specific post process
    virtual void postProcess(
      const std::vector< std::vector<double> >& model_parameter,
      const size_t best_fit_model,
      bool& delete_unused_files) = 0;
    //model-specific tests
    virtual bool testModel(
      const std::vector<double>& parameters) = 0;
  protected:
    GlobalConfig* config;
    SpectralGrid* spectral_grid;
    std::vector<Observation>& observations;

    size_t nb_observation_points = 0;
    size_t nb_spectrum_modifier_param = 0;
    size_t nb_spectral_points = 0;
    
    virtual void setPriors(Priors* priors) = 0;
    virtual void readPriorConfigFile(
      const std::string& file_name, 
      std::vector<std::string>& prior_type, 
      std::vector<std::string>& prior_description, 
      std::vector<std::vector<std::string>>& prior_parameter);

    virtual void convertSpectrumToObservation(
      const std::vector<double>& spectrum, 
      const bool is_flux,
      std::vector<std::vector<double>>& spectrum_obs);

    virtual void convertSpectrumToObservationGPU(
      double* model_spectrum_gpu, 
      const bool is_flux,
      std::vector<double*>& model_spectrum_bands);

    virtual void applyObservationModifier(
      const std::vector<double>& spectrum_modifier_param,
      std::vector<std::vector<double>>& spectrum_obs);
    
    virtual void applyObservationModifierGPU(
      const std::vector<double>& spectrum_modifier_param,
      std::vector<double*>& spectrum_obs);

    virtual void calcPostProcessSpectra(
      const std::vector< std::vector<double> >& model_parameter,
      const size_t best_fit_model,
      std::vector<std::vector< std::vector<double>>>& model_spectra_obs,
      std::vector<double>& spectrum_best_fit);
    virtual void calcPostProcessSpectrum(
      const std::vector<double>& model_parameter,
      std::vector<double>& spectrum,
      std::vector<std::vector<double>>& spectrum_obs);
    virtual void saveBestFitSpectrum(
      const std::vector<double>& spectrum);
    virtual void savePostProcessSpectra(
      const std::vector<std::vector< std::vector<double>>>& model_spectrum_obs);

    virtual bool testCPUvsGPU(
      const std::vector<double>& parameters);
};

}

#endif

