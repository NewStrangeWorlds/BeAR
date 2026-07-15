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
#include "../config/global_config.h"
#include "../spectral_grid/spectral_grid.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../observations/observations.h"
#include "generic_config.h"


namespace bear {



struct ForwardModelOutput{
  bool neglect_model = false;
  std::vector<double> spectrum;
  std::vector<std::vector<double>> spectrum_obs;
};


struct AtmosphereOutput{
  bool neglect_model = false;
  std::vector<double> pressure;
  std::vector<double> altitude;
  std::vector<double> temperature;
  std::vector<std::string> species_symbols;
  std::vector<std::vector<double>> mixing_ratios;
};


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
    virtual bool calcModelCPU(
      const std::vector<double>& physical_parameters, 
      std::vector<double>& spectrum, 
      std::vector<std::vector<double>>& spectrum_obs) = 0;
    //calculate a model on the GPU
    //the return value signals the retrieval to neglect this model
    virtual bool calcModelGPU(
      const std::vector<double>& physical_parameters,
      float* spectrum,
      std::vector<float*>& spectrum_obs) = 0;
    virtual ForwardModelOutput calcModel(
      const std::vector<double>& physical_parameter,
      const bool return_high_res_spectrum);
    virtual AtmosphereOutput getAtmosphereStructure(
      const std::vector<double>& physical_parameter,
      const std::vector<std::string>& species_symbols){
        return AtmosphereOutput();};
    //model-specific post process
    virtual void postProcess(
      const std::vector< std::vector<double> >& posterior_parameters,
      const size_t best_fit_model,
      bool& delete_unused_files) = 0;
    virtual void postProcess(
      GenericConfig* post_process_config_,
      const std::vector< std::vector<double> >& posterior_parameters,
      const size_t best_fit_model,
      bool& delete_unused_files) = 0;
    virtual size_t parametersNumber() = 0;
    //Ordered list of retrieval parameter names for this model's own parameter
    //block (general + module segments), matching the layout of extractParameters.
    //The retrieval layer appends its own tail (high-res Kp/Vsys/dphi/alpha,
    //error inflation) on top of this.
    virtual const std::vector<std::string>& parameterNames() const {
      return parameter_names;}
    virtual void setHighResGrid(SpectralGrid* grid) {
      spectral_grid_highres = grid; }
    //model-specific tests
    virtual bool testModel(
      const std::vector<double>& parameters) = 0;
    const std::vector<double>& spectrumHighRes() const { return spectrum_highres_; }
    float* spectrumHighResGPU() const { return spectrum_highres_gpu_; }
    size_t nbSpectralPointsHighRes() const {
      return spectral_grid_highres ? spectral_grid_highres->nbSpectralPoints() : 0; }
    // Stellar spectrum cached after each calcModel call, used for the per-pixel
    // Doppler correction in high-res kernels (Fs stays at rest; only Fp is shifted).
    // Default returns empty/nullptr for forward models without a stellar spectrum.
    virtual const std::vector<double>& stellarSpectrumCPU() const {
      static const std::vector<double> empty;
      return empty;
    }
    virtual const float* stellarSpectrumGPU() const { return nullptr; }

  protected:
    GlobalConfig* config;
    SpectralGrid* spectral_grid;
    SpectralGrid* spectral_grid_highres = nullptr;
    std::vector<Observation>& observations;

    std::vector<double> spectrum_highres_;
    float* spectrum_highres_gpu_ = nullptr;

    size_t nb_observation_points = 0;
    size_t nb_spectrum_modifier_param = 0;
    size_t nb_spectral_points = 0;

    //Assembled by the derived model (see e.g. EmissionModel::initModules) in the
    //exact order that extractParameters slices the parameter vector.
    std::vector<std::string> parameter_names;

    //Append a sub-module's ordered parameter names to parameter_names. For module
    //types that can be stacked (e.g. multiple cloud layers), pass total>1 to
    //suffix each name with a 1-based instance index so names stay unique across
    //the stack; a single instance keeps its bare name.
    void appendParameterNames(
      const std::vector<std::string>& names, size_t index, size_t total)
    {
      for (const auto& name : names)
        parameter_names.push_back(
          total > 1 ? name + "_" + std::to_string(index + 1) : name);
    }

    virtual void convertSpectrumToObservation(
      const std::vector<double>& spectrum, 
      const bool is_flux,
      std::vector<std::vector<double>>& spectrum_obs);

    virtual void convertSpectrumToObservationGPU(
      float* model_spectrum_gpu,
      const bool is_flux,
      std::vector<float*>& model_spectrum_bands);

    virtual void applyObservationModifier(
      const std::vector<double>& spectrum_modifier_param,
      std::vector<std::vector<double>>& spectrum_obs);
    
    virtual void applyObservationModifierGPU(
      const std::vector<double>& spectrum_modifier_param,
      std::vector<float*>& spectrum_obs);

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

