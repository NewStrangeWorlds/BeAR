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


#include "priors.h"
#include "../observations/observations.h"
#include "../forward_model/forward_model.h"
#include "../forward_model/generic_config.h"

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

  writeParameterLegend();

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

    // Parse the priors file (name defaults to priors.config, overridable via
    // retrieval.toml) into a name-keyed map. The parameter *ordering* is defined
    // by the forward model, not by the file layout, so priors may be listed in
    // any order.
    auto prior_map = Priors::parseConfigToMap(
      config->retrieval_folder_path + config->priors_config_file);

    // Canonical order: the forward model's own parameter block, followed by the
    // retrieval-layer high-res tail (Kp, Vsys, dphi, and optionally alpha).
    std::vector<std::string> ordered_names = forward_model->parameterNames();

    if (has_highres_observations)
    {
      ordered_names.push_back("kp");
      ordered_names.push_back("vsys");
      ordered_names.push_back("dphi");

      // alpha is a free parameter iff the user supplied a prior named "alpha".
      // This replaces the old heuristic of counting lines in priors.config.
      use_free_alpha = prior_map.count("alpha") > 0;

      if (use_free_alpha)
      {
        ordered_names.push_back("alpha");

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

      // Re-injection restores the absolute line amplitude only through a free
      // alpha; with alpha at its maximum-likelihood value the amplitude drops out.
      if (!use_free_alpha)
        for (const auto& obs : highres_observations)
          if (obs.reinjectModel())
            std::cout << "  WARNING: " << obs.observationName()
                      << " uses #reinject_model 1 without an alpha prior; the"
                      << " absolute line amplitude is not constrained\n";
    }

    priors.initFromMap(prior_map, ordered_names);
  }
  catch(std::runtime_error& e)
  {
    std::cout << e.what() << std::endl;
    exit(1);
  }

  setAdditionalPriors();

  priors.printInfo();

  writeParameterLegend();

  if (config->use_gpu)
  {
    for (auto& obs : highres_observations)
      obs.initDeviceMemory();

    initGPUMemory();
  }
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



//Writes the legend that maps posterior columns to prior names next to the other
//retrieval output, once all priors (including any additional ones) are set up.
void Retrieval::writeParameterLegend()
{
  std::string output_folder = config->output_path.empty()
    ? config->retrieval_folder_path : config->output_path;

  if (!output_folder.empty() && output_folder.back() != '/')
    output_folder += '/';

  priors.writeParameterList(output_folder + "cube_parameters.dat");
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
    const auto& stellar_cpu = forward_model->stellarSpectrumCPU();
    const double* stellar_ptr = stellar_cpu.empty() ? nullptr : stellar_cpu.data();

    for (size_t i = 0; i < nb_highres_observations; ++i)
    {
      log_like += highres_observations[i].computeLogLikelihood(
        spectrum_hr, wavelengths_hr, Kp, Vsys, dphi, alpha,
        stellar_ptr, stellar_cpu.size());
    }
  }

  //if the forward model tells us to neglect the current set of parameters,
  //set the likelihood to a low value
  if (neglect == true) log_like = -1e30;


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
      Kp, Vsys, dphi, alpha,
      forward_model->stellarSpectrumGPU());
  }

  //if the forward model tells us to neglect the current set of parameters,
  //set the likelihood to a low value
  if (neglect == true) log_like = -1e30;


  return log_like;
}



std::vector<float> Retrieval::computeHighResModelFlat(
  std::vector<double>& physical_parameters,
  const size_t observation_index,
  const bool apply_projection)
{
  if (!has_highres_observations
      || observation_index >= nb_highres_observations
      || !config->use_gpu)
    return std::vector<float>();

  if (spectral_grid.nbSpectralPoints() > 0)
    initializeOnDevice(spectrum_dev, spectral_grid.nbSpectralPoints());

  for (size_t i=0; i<observations.size(); ++i)
    initializeOnDevice(spectrum_obs_dev[i], observations[i].nbPoints());

  forward_model->calcModelGPU(
    physical_parameters,
    spectrum_dev,
    spectrum_obs_dev);

  const size_t kp_idx = forward_model->parametersNumber();

  return highres_observations[observation_index].modelMatrixGPU(
    forward_model->spectrumHighResGPU(),
    spectral_grid_highres->wavelength_list_gpu,
    forward_model->nbSpectralPointsHighRes(),
    physical_parameters[kp_idx],
    physical_parameters[kp_idx + 1],
    physical_parameters[kp_idx + 2],
    apply_projection,
    forward_model->stellarSpectrumGPU());
}


HighResLayout Retrieval::highResLayout(const size_t observation_index) const
{
  HighResLayout out;

  if (!has_highres_observations || observation_index >= nb_highres_observations)
    return out;

  const auto& obs = highres_observations[observation_index];
  const auto& orders = obs.orders();
  const size_t ne = obs.nbExposures();

  out.valid = true;
  out.nb_orders = obs.nbOrders();
  out.nb_exposures = ne;
  out.kp_ref = obs.kpRef();
  out.vsys_ref = obs.vsysRef();
  out.orbital_phases = obs.orbitalPhases();

  out.order_nb_pixels.reserve(out.nb_orders);
  for (const auto& o : orders)
    out.order_nb_pixels.push_back(static_cast<int>(o.nb_pixels));

  for (const auto& o : orders)
    out.wavelengths.insert(
      out.wavelengths.end(), o.wavelengths.begin(), o.wavelengths.end());

  const auto& ff = obs.filteredFlux();

  if (ff.size() == out.nb_orders)
  {
    out.data_filtered.resize(obs.totalPixels() * ne, 0.0f);
    size_t offset = 0;

    for (size_t ord = 0; ord < out.nb_orders; ++ord)
    {
      const size_t N = orders[ord].nb_pixels;

      for (size_t e = 0; e < ne; ++e)
        for (size_t p = 0; p < N; ++p)
          out.data_filtered[offset * ne + e * N + p] =
            static_cast<float>(ff[ord][e][p]);

      offset += N;
    }
  }

  return out;
}


HighResModelExport Retrieval::computeHighResModelMatrices(
  std::vector<double>& physical_parameters,
  const size_t observation_index,
  const size_t order_index)
{
  HighResModelExport out;

  if (!has_highres_observations
      || observation_index >= nb_highres_observations
      || !config->use_gpu)
    return out;

  const auto& obs = highres_observations[observation_index];

  if (order_index >= obs.nbOrders())
    return out;

  // Build the forward model exactly as logLikelihoodGPU does.
  if (spectral_grid.nbSpectralPoints() > 0)
    initializeOnDevice(spectrum_dev, spectral_grid.nbSpectralPoints());

  for (size_t i=0; i<observations.size(); ++i)
    initializeOnDevice(spectrum_obs_dev[i], observations[i].nbPoints());

  forward_model->calcModelGPU(
    physical_parameters,
    spectrum_dev,
    spectrum_obs_dev);

  const size_t kp_idx = forward_model->parametersNumber();
  const double Kp   = physical_parameters[kp_idx];
  const double Vsys = physical_parameters[kp_idx + 1];
  const double dphi = physical_parameters[kp_idx + 2];

  const std::vector<float> filtered = obs.modelMatrixGPU(
    forward_model->spectrumHighResGPU(),
    spectral_grid_highres->wavelength_list_gpu,
    forward_model->nbSpectralPointsHighRes(),
    Kp, Vsys, dphi, true,
    forward_model->stellarSpectrumGPU());

  const std::vector<float> unfiltered = obs.modelMatrixGPU(
    forward_model->spectrumHighResGPU(),
    spectral_grid_highres->wavelength_list_gpu,
    forward_model->nbSpectralPointsHighRes(),
    Kp, Vsys, dphi, false,
    forward_model->stellarSpectrumGPU());

  const size_t ne = obs.nbExposures();
  const auto& orders = obs.orders();
  const size_t N = orders[order_index].nb_pixels;

  //flat layout is [order_offset*nb_exp + exp*N + pixel]
  size_t offset = 0;
  for (size_t o = 0; o < order_index; ++o)
    offset += orders[o].nb_pixels;

  out.valid = true;
  out.nb_exposures = ne;
  out.nb_pixels = N;
  out.wavelengths = orders[order_index].wavelengths;

  out.model_filtered.assign(ne, std::vector<double>(N, 0.0));
  out.model_unfiltered.assign(ne, std::vector<double>(N, 0.0));

  for (size_t e = 0; e < ne; ++e)
    for (size_t p = 0; p < N; ++p)
    {
      const size_t idx = offset * ne + e * N + p;
      out.model_filtered[e][p]   = filtered[idx];
      out.model_unfiltered[e][p] = unfiltered[idx];
    }

  const auto& ff = obs.filteredFlux();

  if (order_index < ff.size())
    out.data_filtered = ff[order_index];

  return out;
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



std::pair<std::vector<double>, std::vector<double>> Retrieval::convertCubeParameters(
  std::vector<double>& cube)
{
  //convertHypercubeParameters writes the converted values back into the cube,
  //so a too-short cube would corrupt the heap
  if (cube.size() != priors.numberFree())
  {
    std::string error_message =
      "Number of cube parameters ("  + std::to_string(cube.size())
      + ") not equal to the number of free parameters of the forward model ("
      + std::to_string(priors.numberFree()) + ").\n";
    throw InvalidInput(std::string ("Retrieval::convertCubeParameters"), error_message);
  }

  std::vector<double> parameter;
  std::vector<double> physical_parameter;

  convertHypercubeParameters(
    cube.data(),
    cube.size(),
    parameter,
    physical_parameter);

  return std::make_pair(parameter, physical_parameter);
}


void Retrieval::convertHypercubeParameters(
  double *cube,
  const size_t /*nb_free*/,
  std::vector<double>& parameter,
  std::vector<double>& physical_parameter)
{
  const size_t nb_total = priors.number();
  parameter.assign(nb_total, 0.0);
  physical_parameter.assign(nb_total, 0.0);

  for (size_t i = 0; i < nb_total; i++)
  {
    const int j = priors.free_cube_index[i];
    double cube_val;

    if (j < 0)
    {
      cube_val = 0.0;
    }
    else if (priors.distributions[i]->distributionType() == "Linked prior")
    {
      const int src_j = priors.free_cube_index[priors.prior_links[i]];
      cube_val = (src_j >= 0) ? cube[src_j] : 0.0;
    }
    else
    {
      cube_val = cube[j];
    }

    parameter[i]          = priors.distributions[i]->parameterValue(cube_val);
    physical_parameter[i] = priors.distributions[i]->parameterPhysicalValue(cube_val);
  }

  for (size_t i = 0; i < nb_total; i++)
  {
    const int j = priors.free_cube_index[i];
    if (j >= 0) cube[j] = parameter[i];
  }
}


std::vector<double> Retrieval::convertToPhysicalParameters(
  const std::vector<double>& parameters)
{
  if (parameters.size() != priors.numberFree())
  {
    std::string error_message =
      "Number of posterior parameters not equal to the number of free parameters of the forward model.\n";
    throw InvalidInput(std::string ("Retrieval::convertToPhysicalParameters"), error_message);
  }

  const size_t nb_total = priors.number();
  std::vector<double> physical_parameters(nb_total, 0.0);

  for (size_t i = 0; i < nb_total; ++i)
  {
    const int j = priors.free_cube_index[i];
    if (j < 0)
      physical_parameters[i] = priors.distributions[i]->parameterPhysicalValue(0.0);
    else
      physical_parameters[i] = priors.distributions[i]->applyParameterUnit(parameters[j]);
  }

  return physical_parameters;
}


std::vector<double> Retrieval::convertToParameterUnits(
  const std::vector<double>& physical_parameters)
{
  if (physical_parameters.size() != priors.numberFree())
  {
    std::string error_message =
      "Number of physical parameters not equal to the number of free parameters of the forward model.\n";
    throw InvalidInput(std::string ("Retrieval::convertToParameterUnits"), error_message);
  }

  const size_t nb_total = priors.number();
  std::vector<double> parameters(priors.numberFree(), 0.0);

  for (size_t i = 0; i < nb_total; ++i)
  {
    const int j = priors.free_cube_index[i];
    if (j >= 0 && priors.distributions[i]->distributionType() != "Linked prior")
      parameters[j] = priors.distributions[i]->invertParameterUnit(physical_parameters[j]);
  }

  return parameters;
}


}
