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


#include "retrieval.h"


#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <csignal>
#include <cstdlib>


#include "multinest_parameter.h"
#include "priors.h"
#include "../observations/observations.h"
#include "../forward_model/forward_model.h"
#include "../forward_model/generic_config.h"

#include "../../_deps/multinest-src/MultiNest_v3.12_CMake/multinest/include/multinest.h"
#include "../CUDA_kernels/data_management_kernels.h"
#include "../additional/exceptions.h"

#include "../forward_model/transmission/transmission.h"


namespace bear{


bool stop_model = false;

void signalHandler(int sig) 
{
  std::cout << "Received signal " << sig << "\n"; 
  
  if (sig == SIGCONT) 
    stop_model = true;

}


Retrieval::Retrieval(GlobalConfig* global_config) 
  : Retrieval(global_config, std::string(""))
{
  
}



Retrieval::Retrieval(
  GlobalConfig* global_config,
  GenericConfig* model_config,
  const std::vector<ObservationInput>& observation_input,
  const std::vector<PriorConfig>& prior_config)
  : spectral_grid(global_config)
{
  config = global_config;

  size_t nb_add_priors = 0;

  if (config->use_error_inflation)
    nb_add_priors += 1;

  try
  {
    setObservations(observation_input);
    
    std::cout << "\nTotal number of wavelength points: " 
              << spectral_grid.nbSpectralPoints() << "\n\n";
    
    forward_model = selectForwardModel(config->forward_model_type, model_config);
    
    priors.init(
      prior_config, 
      forward_model->parametersNumber() + nb_add_priors);
  }
  catch(std::runtime_error& e) 
  {
    std::cout << e.what() << std::endl;
    exit(1);
  }

  priors.printInfo();

  if (config->use_gpu)
    initGPUMemory();
}



Retrieval::Retrieval(
  GlobalConfig* global_config,
  const std::string additional_observation_file) 
  : spectral_grid(global_config)
{ 
  config = global_config;

  std::signal(SIGCONT, signalHandler);

  std::string folder = config->retrieval_folder_path;
  std::string observation_folder = folder;

  //try to initialise the model
  //if there is an error, we exit the retrieval
  try
  {
    std::vector<std::string> file_list, modifier_list;

    loadObservationFileList(
      observation_folder, 
      file_list,
      modifier_list);

    //if we do postprocessing, we may need to read in the file that describes the maximum wavelength range
    //spectra will be generated for
    //this is necessary to obtain an estimate for the effective temperature
    if (additional_observation_file.size() > 0)
    {
      std::string postprocess_spectrum_data = config->retrieval_folder_path + additional_observation_file;
      std::fstream file(postprocess_spectrum_data.c_str(), std::ios::in);

      if (!file.fail())
      {
        file_list.push_back(additional_observation_file);
        modifier_list.push_back("none");
        file.close(); 
      }
    }

    if (!file_list.empty())
    {
      loadObservations(
        observation_folder,
        file_list,
        modifier_list);

      std::cout << "\nTotal number of low-res wavelength points: "
                << spectral_grid.nbSpectralPoints() << "\n\n";
    }

    // Load high-res observations and create the high-res spectral grid
    loadHighResObservations(observation_folder);
    
    forward_model = selectForwardModel(config->forward_model_type, nullptr);
    
    // Set up dual spectral grid if high-res observations exist
    if (has_highres_observations)
      forward_model->setHighResGrid(spectral_grid_highres.get());

    // High-res parameters: always Kp and Vsys; optionally alpha.
    // Detect alpha by counting lines in priors.config vs expected model params.
    size_t nb_highres_param = 0;

    if (has_highres_observations)
    {
      nb_highres_param = 3;  // Kp, Vsys, dphi

      // Peek at priors.config to check if alpha is included (3rd high-res param).
      // Count lines the same way as Priors::readConfigFile: every non-empty line.
      const std::string priors_file = config->retrieval_folder_path + "priors.config";
      std::ifstream pf(priors_file);
      size_t nb_prior_lines = 0;
      std::string line;

      while (std::getline(pf, line))
      {
        if (!line.empty())
          ++nb_prior_lines;
      }

      if (nb_prior_lines == forward_model->parametersNumber() + 4)
      {
        use_free_alpha = true;
        nb_highres_param = 4;

        for (auto& obs : highres_observations)
        {
          if (obs.hasFluxUncertainties())
          {
            obs.likelihood_mode = HighResLikelihoodMode::gibson;
            obs.precomputeGibsonStatistics();
            std::cout << "  Using Gibson Eq. 4 likelihood (per-pixel uncertainties)\n";
          }
          else
          {
            obs.likelihood_mode = HighResLikelihoodMode::free_alpha;
          }
        }

        std::cout << "  Alpha is a free retrieval parameter\n";
      }
    }

    priors.init(
      config->retrieval_folder_path,
      forward_model->parametersNumber() + nb_highres_param);
  }
  catch(std::runtime_error& e)
  {
    std::cout << e.what() << std::endl;
    exit(1);
  }

  setAdditionalPriors();

  priors.printInfo();

  if (config->use_gpu)
  {
    for (auto& obs : highres_observations)
      obs.initDeviceMemory();

    initGPUMemory();
  }
}




bool Retrieval::run()
{
  //Configure Multinest
  MultinestParameter param(config);

  size_t nb_free = priors.numberFree();

  param.ndims = nb_free;
  param.nPar = nb_free;
  param.nClsPar = nb_free;

  for (size_t i = 0; i < nb_free; ++i)
    param.pWrap[i] = 0;

  //We give the MultiNest function a pointer to the retrieval class
  //That way, we can access the current retrieval object and its static member routines
  param.context = this;


  //Call MultiNest
  if (config->use_gpu == false)
    nested::run(
      param.is,
      param.mmodal,
      param.ceff,
      param.nlive,
      param.tol,
      param.efr,
      param.ndims,
      param.nPar,
      param.nClsPar,
      param.maxModes,
      param.updInt,
      param.Ztol,
      param.root,
      param.seed,
      param.pWrap,
      param.fb,
      param.resume,
      param.outfile,
      param.initMPI,
      param.logZero,
      param.maxiter,
      Retrieval::multinestLogLike,
      Retrieval::multinestDumper,
      param.context);
  else
    nested::run(
      param.is,
      param.mmodal,
      param.ceff,
      param.nlive,
      param.tol,
      param.efr,
      param.ndims,
      param.nPar,
      param.nClsPar,
      param.maxModes,
      param.updInt,
      param.Ztol,
      param.root,
      param.seed,
      param.pWrap,
      param.fb,
      param.resume,
      param.outfile,
      param.initMPI,
      param.logZero,
      param.maxiter,
      Retrieval::multinestLogLikeGPU,
      Retrieval::multinestDumper,
      param.context);

  return true;
}



void Retrieval::setAdditionalPriors()
{
  if (config->use_error_inflation)
  {
    //this creates the prior distribution for the error exponent
    //first, we need to find the minimum and maximum values of the observational data errors
    double error_max = 0;

    for (auto & obs : observations)
    {
      double obs_error_max = *std::max_element(
        std::begin(obs.data_error),
        std::end(obs.data_error));

      if (obs_error_max > error_max)
        error_max = obs_error_max;
    }

    double error_min = error_max;

    for (auto & obs : observations)
    {
      double obs_error_min = *std::min_element(
        std::begin(obs.data_error),
        std::end(obs.data_error));

      if (obs_error_min < error_min)
        error_min = obs_error_min;
    }

    error_min = std::log10(0.1 * error_min * error_min);
    error_max = std::log10(100.0 * error_max * error_max);

    priors.add(
      std::vector<PriorConfig> {
        PriorConfig(
          std::string("uniform"),
          std::string("error exponent"),
          std::vector<double>{error_min, error_max})});
  }
}



double Retrieval::computeLikelihood(
  std::vector<double>& physical_parameter)
{
  double log_like = 0;

  if (!config->use_gpu)
  {
    log_like =  logLikelihood(physical_parameter);
  }
  else
  {
    log_like =  logLikelihoodGPU(physical_parameter);
  }

  return log_like;
}



double Retrieval::logLikelihood(
  std::vector<double>& physical_parameters)
{
  std::vector<double> model_spectrum(
    spectral_grid.nbSpectralPoints(),
    0.0);

  std::vector<std::vector<double>> model_spectrum_obs(
    nb_observations,
    std::vector<double>{});


  bool neglect = forward_model->calcModelCPU(
    physical_parameters,
    model_spectrum,
    model_spectrum_obs);


  double error_inflation = 0;

  if (config->use_error_inflation)
    error_inflation = std::pow(10, physical_parameters.back());


  double log_like = 0;

  // Low-res chi-square likelihood
  for (size_t i=0; i<observations.size(); ++i)
  {
    for (size_t j=0; j<observations[i].nbPoints(); ++j)
    {
      //Eq. 22 from Paper I
      const double error_square =
        observations[i].data_error[j]
        * observations[i].data_error[j]
        + error_inflation;

      const double obs_delta = observations[i].data[j] - model_spectrum_obs[i][j];

      //Eq. 23 from Paper I
      log_like +=
        (- 0.5 * std::log(error_square* 2.0 * constants::pi)
         - 0.5 * obs_delta*obs_delta / error_square)
         * observations[i].likelihood_weight[j];
    }
  }

  // High-res Brogi & Line 2019 likelihood
  if (has_highres_observations)
  {
    const size_t kp_idx = forward_model->parametersNumber();
    const double Kp   = physical_parameters[kp_idx];
    const double Vsys = physical_parameters[kp_idx + 1];
    const double dphi = physical_parameters[kp_idx + 2];
    const double alpha = use_free_alpha ? physical_parameters[kp_idx + 3] : 1.0;

    const auto& spectrum_hr = forward_model->spectrumHighRes();
    const auto& wavelengths_hr = spectral_grid_highres->wavelength_list;

    for (size_t i = 0; i < nb_highres_observations; ++i)
    {
      log_like += highres_observations[i].computeLogLikelihood(
        spectrum_hr, wavelengths_hr, Kp, Vsys, dphi, alpha);
    }
  }

  //if the forward model tells us to neglect the current set of parameters,
  //set the likelihood to a low value
  if (neglect == true) log_like = -1e30;


  if (config->multinest_print_iter_values)
    std::cout << log_like << "\n";

  return log_like;
}




double Retrieval::logLikelihoodGPU(
  std::vector<double>& physical_parameters)
{
  if (spectral_grid.nbSpectralPoints() > 0)
    initializeOnDevice(
      spectrum_dev,
      spectral_grid.nbSpectralPoints());

  for (size_t i=0; i<observations.size(); ++i)
    initializeOnDevice(
      spectrum_obs_dev[i],
      observations[i].nbPoints());


  bool neglect = forward_model->calcModelGPU(
    physical_parameters,
    spectrum_dev,
    spectrum_obs_dev);


  double error_inflation = 0;

  if (config->use_error_inflation)
    error_inflation = std::pow(10, physical_parameters.back());


  double log_like = logLikeDev(spectrum_obs_dev, error_inflation);

  // High-res Brogi & Line 2019 likelihood (fully on GPU)
  if (has_highres_observations)
  {
    const size_t kp_idx = forward_model->parametersNumber();
    const double Kp   = physical_parameters[kp_idx];
    const double Vsys = physical_parameters[kp_idx + 1];
    const double dphi = physical_parameters[kp_idx + 2];
    const double alpha = use_free_alpha ? physical_parameters[kp_idx + 3] : 1.0;

    log_like += logLikeHighResDev(
      forward_model->spectrumHighResGPU(),
      spectral_grid_highres->wavelength_list_gpu,
      forward_model->nbSpectralPointsHighRes(),
      Kp, Vsys, dphi, alpha);
  }

  //if the forward model tells us to neglect the current set of parameters,
  //set the likelihood to a low value
  if (neglect == true) log_like = -1e30;


  if (config->multinest_print_iter_values)
    std::cout << log_like << "\n";

  return log_like;
}



ForwardModelOutput Retrieval::computeModel(
  std::vector<double>& physical_parameters,
  const bool return_high_res_spectrum)
{
  return forward_model->calcModel(
    physical_parameters, 
    return_high_res_spectrum);
}


AtmosphereOutput Retrieval::computeAtmosphereStructure(
  std::vector<double>& physical_parameters,
  const std::vector<std::string>& species_symbols)
{
  return forward_model->getAtmosphereStructure(
    physical_parameters, 
    species_symbols);
}



void Retrieval::initGPUMemory()
{
  if (gpu_memory_initialized)
    return;

  if (spectral_grid.nbSpectralPoints() > 0)
    allocateOnDevice(spectrum_dev, spectral_grid.nbSpectralPoints());

  spectrum_obs_dev.resize(observations.size(), nullptr);

  for (size_t i=0; i<observations.size(); ++i)
    allocateOnDevice(spectrum_obs_dev[i], observations[i].nbPoints());

  allocateOnDevice(d_log_like_dev, 1);

  gpu_memory_initialized = true;
}


void Retrieval::freeGPUMemory()
{
  if (!gpu_memory_initialized)
    return;

  if (spectrum_dev != nullptr)
    deleteFromDevice(spectrum_dev);

  for (size_t i=0; i<spectrum_obs_dev.size(); ++i)
    deleteFromDevice(spectrum_obs_dev[i]);

  spectrum_obs_dev.clear();

  deleteFromDevice(d_log_like_dev);

  gpu_memory_initialized = false;
}


Retrieval::~Retrieval()
{
  freeGPUMemory();
}


}
