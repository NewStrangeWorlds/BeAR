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
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "phase_curve.h"

#include "../../chemistry/chem_species.h"
#include "../../CUDA_kernels/data_management_kernels.h"
#include "../../CUDA_kernels/contribution_function_kernels.h"
#include "../atmosphere/atmosphere.h"
#include "../../additional/aux_functions.h"


namespace bear{


PhaseCurvePostProcessConfig::PhaseCurvePostProcessConfig (
  const std::string& folder_path,
  const std::string& file_name)
{
  readConfigFile(folder_path + file_name);
}



PhaseCurvePostProcessConfig::PhaseCurvePostProcessConfig (
  const bool save_temperatures_,
  const bool save_spectra_,
  const bool save_contribution_functions_,
  const std::vector<std::string>& species_to_save_)
{
  save_temperatures = save_temperatures_;
  save_spectra = save_spectra_;
  save_contribution_functions = save_contribution_functions_;

  species_to_save = aux::findChemicalSpecies(species_to_save_);
}



void PhaseCurvePostProcessConfig::readConfigFile(const std::string& file_name)
{
  std::ifstream test_file(file_name.c_str());

  if (!test_file.good())
  {
    std::cout << "\n Post-process config file not found. Using default options!\n\n";

    return;
  }
  test_file.close();

  std::cout << "\nParameters read from " << file_name << " :\n";

  toml::table cfg = parseConfigFile(file_name);

  delete_sampler_files = readBooleanParameter(cfg, "delete_sampler_files", delete_sampler_files);
  save_spectra = readBooleanParameter(cfg, "save_spectra", save_spectra);
  save_temperatures = readBooleanParameter(cfg, "save_temperatures", save_temperatures);
  save_contribution_functions = readBooleanParameter(cfg, "save_contribution_functions", save_contribution_functions);
  species_to_save = readChemicalSpecies(cfg, "species_to_save");
  //BeAR-highres addition: dump the high-res spectrum for every posterior
  //sample (spectra_hr_all.bin), needed for the per-species decomposition figure
  save_hr_spectra_all = readBooleanParameter(cfg, "save_hr_spectra_all", save_hr_spectra_all);
}


void PhaseCurveModel::postProcess(
  GenericConfig* post_process_config_,
  const std::vector< std::vector<double> >& model_parameter,
  const size_t best_fit_model,
  bool& delete_unused_files)
{
  PhaseCurvePostProcessConfig post_process_config =
    *(dynamic_cast<PhaseCurvePostProcessConfig*>(post_process_config_));

  if (post_process_config.delete_sampler_files)
    delete_unused_files = true;

  postProcess(
    post_process_config,
    model_parameter,
    best_fit_model);
}



void PhaseCurveModel::postProcess(
  const std::vector< std::vector<double> >& model_parameter,
  const size_t best_fit_model,
  bool& delete_unused_files)
{
  PhaseCurvePostProcessConfig post_process_config(
    config->retrieval_folder_path, config->post_process_config_file);

  if (post_process_config.delete_sampler_files)
    delete_unused_files = true;

  postProcess(
    post_process_config,
    model_parameter,
    best_fit_model);
}


void PhaseCurveModel::postProcess(
  const PhaseCurvePostProcessConfig& post_process_config,
  const std::vector< std::vector<double> >& model_parameter,
  const size_t best_fit_model)
{
  const size_t nb_models = model_parameter.size();

  if (post_process_config.save_spectra)
  {
    std::vector<std::vector<std::vector<double>>> model_spectra_obs;

    std::vector<double> model_spectrum_best_fit;

    save_hr_spectra_all_ = post_process_config.save_hr_spectra_all;

    calcPostProcessSpectra(
      model_parameter,
      best_fit_model,
      model_spectra_obs,
      model_spectrum_best_fit);

    saveBestFitSpectrum(model_spectrum_best_fit);
    savePostProcessSpectra(model_spectra_obs);
  }

  std::vector<std::vector<double>> temperature_profiles(
    nb_models,
    std::vector<double>(nb_grid_points, 0));

  std::vector<std::vector<std::vector<double>>> mixing_ratios(
    nb_models,
    std::vector<std::vector<double>>(constants::species_data.size(),
    std::vector<double>(nb_grid_points,0)));

  for (size_t i=0; i<nb_models; ++i)
  {
    postProcessModel(
      model_parameter[i],
      temperature_profiles[i],
      mixing_ratios[i]);

    if (i == best_fit_model && post_process_config.save_contribution_functions)
      postProcessContributionFunctions(model_parameter[i]);
  }

  if (post_process_config.species_to_save.size() > 0)
    for (auto & i : post_process_config.species_to_save)
      savePostProcessChemistry(mixing_ratios, i);

  if (post_process_config.save_temperatures)
    savePostProcessTemperatures(temperature_profiles);
}



void PhaseCurveModel::postProcessModel(
  const std::vector<double>& parameters,
  std::vector<double>& temperature_profile,
  std::vector<std::vector<double>>& mixing_ratios)
{
  extractParameters(parameters);

  calcAtmosphereStructure(parameters);

  for (auto & i : constants::species_data)
  {
    for (size_t j=0; j<nb_grid_points; ++j)
      mixing_ratios[i.id][j] = atmosphere.number_densities[j][i.id]/atmosphere.number_densities[j][_TOTAL];
  }

  temperature_profile = atmosphere.temperature;
}



void PhaseCurveModel::savePostProcessChemistry(
  const std::vector<std::vector<std::vector<double>>>& mixing_ratios,
  const unsigned int species)
{
  std::fstream file;
  std::string file_name = config->post_output_path + "/chem_";

  file_name += constants::species_data[species].symbol;
  file_name += ".dat";

  file.open(file_name.c_str(), std::ios::out);


  const size_t nb_models = mixing_ratios.size();

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    file << std::setprecision(10) << std::scientific << atmosphere.pressure[i];

    for (size_t j=0; j<nb_models; ++j)
      file << "\t" << mixing_ratios[j][species][i];

    file << "\n";
  }

}



void PhaseCurveModel::savePostProcessTemperatures(
  const std::vector<std::vector<double>>& temperature_profiles)
{
  std::fstream file;
  std::string file_name = config->post_output_path + "/temperature_structures.dat";
  file.open(file_name.c_str(), std::ios::out);

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    file << std::setprecision(10) << std::scientific << atmosphere.pressure[i];

    for(size_t j=0; j<temperature_profiles.size(); ++j)
      file << "\t" << temperature_profiles[j][i];

    file << "\n";
  }

}


void PhaseCurveModel::postProcessContributionFunctions(
  const std::vector<double>& parameter)
{
  //Contribution functions are computed on the LOW-RES grid and written per
  //low-res observation. A high-res-only run (e.g. IGRINS alone) has neither,
  //so the opacity kernel below would be launched with zero spectral points
  //and abort with "invalid argument" -- at the very end of post-processing,
  //after all the useful output has already been produced.
  //To get contribution functions for a high-res-only retrieval, add a
  //postprocess_spectrum_data.dat to the run folder: Retrieval appends it as a
  //low-res observation, which gives this path a grid to work on.
  if (observations.empty() || spectral_grid->nbSpectralPoints() == 0)
  {
    std::cout << "\nSkipping contribution functions: this retrieval has no "
              << "low-resolution observations.\n"
              << "Add a postprocess_spectrum_data.dat to compute them.\n";
    return;
  }

  std::vector<double> cloud_parameters(
      parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param,
      parameter.begin() + nb_general_param + nb_total_chemistry_param + nb_temperature_param + nb_total_cloud_param);

  opacity_calc.calculateGPU(cloud_models, cloud_parameters);

  float* contribution_functions_dev = nullptr;
  size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  allocateOnDevice(contribution_functions_dev, nb_spectral_points*nb_grid_points);

  contributionFunctionGPU(
    contribution_functions_dev,
    opacity_calc.absorption_coeff_gpu,
    spectral_grid->wavenumber_list_gpu,
    atmosphere.temperature,
    atmosphere.altitude,
    nb_spectral_points);

  std::vector<float> contribution_functions_float(nb_spectral_points*nb_grid_points, 0.0f);
  moveToHost(contribution_functions_dev, contribution_functions_float);
  deleteFromDevice(contribution_functions_dev);

  std::vector<double> contribution_functions_all(contribution_functions_float.begin(), contribution_functions_float.end());

  std::vector< std::vector<double> > contribution_functions(nb_grid_points, std::vector<double>(nb_spectral_points, 0));

  for (size_t i=0; i<nb_spectral_points; ++i)
    for (size_t j=0; j<nb_grid_points; ++j)
      contribution_functions[j][i] = contribution_functions_all[j*nb_spectral_points + i];

  for (size_t i=0; i<observations.size(); ++i)
  {
    std::vector<std::vector<double>> contribution_functions_band(
      nb_grid_points,
      std::vector<double>(observations[i].nbPoints(), 0));

    const bool is_flux = false;

    for (size_t j=0; j<nb_grid_points; ++j)
      contribution_functions_band[j] =
        observations[i].processModelSpectrum(contribution_functions[j], is_flux);

    saveContributionFunctions(contribution_functions_band, i);
  }

}



void PhaseCurveModel::saveContributionFunctions(
  std::vector< std::vector<double>>& contribution_function, const size_t observation_index)
{
  std::string observation_name = observations[observation_index].observationName();
  std::replace(observation_name.begin(), observation_name.end(), ' ', '_');

  std::string file_name = config->post_output_path + "/contribution_function_" + observation_name + ".dat";


  std::fstream file(file_name.c_str(), std::ios::out);

  for (size_t j=0; j<nb_grid_points; ++j)
  {
    file << std::setprecision(10) << std::scientific << atmosphere.pressure[j] << "\t";

    for (size_t i=0; i<observations[observation_index].spectral_bands.nbBands(); ++i)
      file << std::setprecision(10) << std::scientific << contribution_function[j][i] << "\t";

    file << "\n";
  }

  file.close();
}


// Override the base-class calcPostProcessSpectra to read back spectrum_highres_gpu_
// for the best-fit sample only. For HRCCS retrievals there is no meaningful per-sample
// high-res output: the data lives in CCF/Kp-Vsys space, not in flux space.
void PhaseCurveModel::calcPostProcessSpectra(
  const std::vector<std::vector<double>>& model_parameter,
  const size_t best_fit_model,
  std::vector<std::vector<std::vector<double>>>& model_spectra_obs,
  std::vector<double>& spectrum_best_fit)
{
  const size_t nb_models = model_parameter.size();
  const size_t nb_hr = nbSpectralPointsHighRes();
  const bool has_hr = (nb_hr > 0 && spectrum_highres_gpu_ != nullptr);
  const bool dump_hr = (save_hr_spectra_all_ && has_hr);

  model_spectra_obs.resize(nb_models);

  std::ofstream hr_bin_file;
  if (dump_hr)
  {
    hr_bin_file.open(
      config->post_output_path + "/spectra_hr_all.bin",
      std::ios::binary);
    const int32_t n_m = static_cast<int32_t>(nb_models);
    const int32_t n_w = static_cast<int32_t>(nb_hr);
    hr_bin_file.write(reinterpret_cast<const char*>(&n_m), sizeof(int32_t));
    hr_bin_file.write(reinterpret_cast<const char*>(&n_w), sizeof(int32_t));

    std::ofstream wl_file(config->post_output_path + "/spectra_hr_wavelengths.dat");
    for (size_t k = 0; k < nb_hr; ++k)
      wl_file << std::setprecision(10) << std::scientific
              << spectral_grid_highres->wavelength_list[k] << "\n";
  }

  std::cout << "\n";

  for (size_t i = 0; i < nb_models; ++i)
  {
    std::cout << "\rPostprocess spectra, model " << i << " of " << nb_models << std::flush;

    std::vector<double> spectrum_lowres;
    calcPostProcessSpectrum(model_parameter[i], spectrum_lowres, model_spectra_obs[i]);

    if (dump_hr)
    {
      std::vector<float> hr_float(nb_hr);
      moveToHost(spectrum_highres_gpu_, hr_float);
      hr_bin_file.write(reinterpret_cast<const char*>(hr_float.data()),
                        static_cast<std::streamsize>(nb_hr * sizeof(float)));
    }

    if (i == best_fit_model)
    {
      if (has_hr)
      {
        std::vector<float> hr_float(nb_hr);
        moveToHost(spectrum_highres_gpu_, hr_float);
        spectrum_best_fit.assign(hr_float.begin(), hr_float.end());
      }
      else
      {
        spectrum_best_fit = spectrum_lowres;
      }
    }
  }

  std::cout << "\n";
}


// Override saveBestFitSpectrum: if the spectrum matches the high-res grid size,
// use the high-res wavelength list instead of the low-res one.
void PhaseCurveModel::saveBestFitSpectrum(const std::vector<double>& spectrum)
{
  if (spectral_grid_highres &&
      spectrum.size() == spectral_grid_highres->nbSpectralPoints())
  {
    const std::string file_name = config->post_output_path + "/spectrum_best_fit_hr.dat";
    std::fstream file(file_name.c_str(), std::ios::out);
    for (size_t i = 0; i < spectrum.size(); ++i)
      file << std::setprecision(10) << std::scientific
           << spectral_grid_highres->wavelength_list[i] << "\t"
           << spectrum[i] << "\n";
    return;
  }
  ForwardModel::saveBestFitSpectrum(spectrum);
}


}
