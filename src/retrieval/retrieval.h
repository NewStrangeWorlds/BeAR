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


#ifndef _retrieval_h
#define _retrieval_h

#include <vector>
#include <string>
#include <iostream>
#include <memory>


#include "../spectral_grid/spectral_grid.h"
#include "../observations/observations.h"
#include "../observations/highres_observation.h"
#include "../config/global_config.h"
#include "priors.h"


namespace bear {

extern bool stop_model;

void signalHandler(int sig);

//forward declaration
class ForwardModel;
struct ForwardModelOutput;
struct AtmosphereOutput;
class GenericConfig;


//the main class that does the retrieval
class Retrieval{
  public:
    Retrieval(
      GlobalConfig* global_config, 
      const std::string additional_observation_file);
    Retrieval(GlobalConfig* global_config);
    Retrieval(
      GlobalConfig* global_config,
      GenericConfig* model_config,
      const std::vector<ObservationInput>& observation_input,
      const std::vector<PriorConfig>& prior_config);
    ~Retrieval();
    
    GlobalConfig* config = nullptr;
    SpectralGrid spectral_grid;
    std::unique_ptr<SpectralGrid> spectral_grid_highres;
    std::vector<Observation> observations;
    std::vector<HighResObservation> highres_observations;
    Priors priors;

    size_t nb_observations = 0;
    size_t nb_highres_observations = 0;
    bool has_highres_observations = false;
    bool use_free_alpha = false;  // when true, alpha is a free retrieval parameter

    virtual bool run();

    std::pair<std::vector<double>, std::vector<double>> convertCubeParameters(
      std::vector<double>& cube);
    std::vector<double> convertToPhysicalParameters(
      const std::vector<double>& parameters);
    
    double computeLikelihood(
      std::vector<double>& parameters);

    ForwardModelOutput computeModel(
      std::vector<double>& physical_parameters,
      const bool return_high_res_spectrum);
    
    AtmosphereOutput computeAtmosphereStructure(
      std::vector<double>& physical_parameters,
      const std::vector<std::string>& species_symbols);
    
    size_t nbParameters() {
      return priors.numberFree();}

  protected:
    std::unique_ptr<ForwardModel> forward_model;

    std::unique_ptr<ForwardModel> selectForwardModel(
       const std::string model_description,
       GenericConfig* model_config);
    void setAdditionalPriors();
    
    void setObservations(
      const std::vector<ObservationInput>& observation_input);
    void loadObservations(
      const std::string file_folder, 
      const std::vector<std::string>& file_list,
      const std::vector<std::string>& modifier_list);
    void loadObservationFileList(
      const std::string file_folder,
      std::vector<std::string>& file_list,
      std::vector<std::string>& modifier_list);
    void loadHighResObservations(const std::string& file_folder);
    std::vector<std::string> highres_file_list;
  private:
    double logLikelihood(
      std::vector<double>& parameters);
    double logLikelihoodGPU(
      std::vector<double>& parameters);
    void convertHypercubeParameters(
      double *cube,
      const size_t nb_param,
      std::vector<double>& parameter,
      std::vector<double>& physical_parameter);
    static void multinestLogLike(
      double *Cube, 
      int &ndim, 
      int &npars, 
      double &lnew, 
      void *context);
    static void multinestLogLikeGPU(
      double *cube, 
      int &nb_dim, 
      int &nb_param, 
      double &new_log_like, 
      void *context);
    double logLikeDev(
      std::vector<float*> model_spectrum,
      const double error_inflation_coefficient);
    double logLikeHighResDev(
      float* spectrum_hr_gpu,
      double* model_wl_gpu,
      size_t nb_hr_points,
      double Kp, double Vsys, double dphi, double alpha);

    float* spectrum_dev = nullptr;
    std::vector<float*> spectrum_obs_dev;
    double* d_log_like_dev = nullptr;
    bool gpu_memory_initialized = false;

    void initGPUMemory();
    void freeGPUMemory();

    static void multinestDumper(
      int &nSamples, 
      int &nlive, 
      int &nPar, 
      double **physLive, 
      double **posterior,
      double **paramConstr, 
      double &maxLogLike, 
      double &logZ, 
      double &INSlogZ, 
      double &logZerr, 
      void *context);
};


}


#endif
