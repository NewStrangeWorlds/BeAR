import os
import sys

current_directory = os.path.dirname(os.path.realpath(__file__))
parent_directory = os.path.dirname(current_directory)
sys.path.append(parent_directory)

from lib import pybear
import numpy as np
import pymultinest


#setting the basic properties of the model
use_gpu = True
model_type = "transmission"
cross_section_file_path = "/media/data/opacity_data/helios-k/"
spectral_discretisation = 'const_wavenumber'
resolution = 1.0

retrieval_folder = "../../TransmissionExample/"
multinest_output_folder = "../../TransmissionExample/"
post_output_folder = "../../TransmissionExample/"

#create the general model config
model_config = pybear.Config(
  use_gpu, 
  model_type, 
  cross_section_file_path, 
  spectral_discretisation, 
  resolution,
  multinest_output_folder,
  post_output_folder)

#read in the observation data
observation = np.loadtxt(retrieval_folder + "WASP-12b_kreidberg.dat", skiprows=11)
obs_wavelength_bins = observation[:,0:2]
obs_data = observation[:,2]
obs_error = observation[:,3]

#and create a pyBeAR observation input
wasp_12b_obs = pybear.Observation(
  "WASP-12b_Kreidberg", "band-spectroscopy", obs_wavelength_bins, obs_data, obs_error)

#create a list of all observations
observations = list([wasp_12b_obs])


#create the list of priors
priors_config = list([
  pybear.Prior("gaussian", "log_g", [2.99, 0.03]),
  pybear.Prior("uniform", "r_planet", [1.6, 2.3], "Rj"),
  pybear.Prior("gaussian", "r_star", [1.57, 0.07], "Rs"),
  pybear.Prior("log_uniform", "H2O", [1e-12, 1e-1]),
  pybear.Prior("log_uniform", "TiO", [1e-12, 1e-3]),
  pybear.Prior("log_uniform", "VO", [1e-12, 1e-3]),
  pybear.Prior("log_uniform", "K", [1e-12, 1e-3]),
  pybear.Prior("uniform", "temperature", [500, 1600])])


#Now we configure the forward model
nb_grid_points = 100
bottom_pressure = 10
top_pressure = 1e-5

opacity_species_data = np.array([
  ['CIA-H2-H2', 'CIA/H2-H2'], 
  ['CIA-H2-He', 'CIA/H2-He'],
  ['H2', 'Rayleigh'],
  ['He', 'Rayleigh'],
  ['H2O', 'Molecules/1H2-16O__POKAZATEL_e2b_n'],
  ['CO', 'Molecules/12C-16O__Li2015_e2b_n'],
  ['CO2', 'Molecules/12C-16O2__UCL-4000_e2b'],
  ['H2O', 'Rayleigh'],
  ['CO', 'Rayleigh'],
  ['CO2', 'Rayleigh']])
  
#configuration of the temperature profile
temperature_profile = "const"
#a constant temperature profile does not need any config parameters
temperature_profile_parameters = [] 

#configuration of the chemistry
#here, we use two: one with isoprofiles and the background chemistry
chemistry_models = [
  "iso", 
  "bg"]
#the chemical species for each of the two models
chemistry_parameters = [
  ["H2O", "TiO", "VO", "K"], 
  ["H2He"]]

opacity_species_symbols = opacity_species_data[:, 0]
opacity_species_folders = opacity_species_data[:, 1]

#Create the configuration of the forward model
forward_model_config = pybear.TransmissionModelConfig(
  nb_grid_points,
  bottom_pressure,
  top_pressure,
  temperature_profile,
  temperature_profile_parameters,
  chemistry_models,
  chemistry_parameters,
  opacity_species_symbols,
  opacity_species_folders)


#Now, we can create the BeAR retrieval object
model = pybear.Retrieval(
  model_config, 
  forward_model_config, 
  observations,
  priors_config)


#Define the priors and likelihood functions for MultiNest
#We use the internal functions of BeAR to convert the cube parameters
#and to compute the likelihood
cube_parameters = np.zeros(model.nbParameters())
parameters = np.zeros(model.nbParameters())
physical_parameters = np.zeros(model.nbParameters())


def priors(cube, ndim, nparams):
  for i in range(nparams):
    cube_parameters[i] = cube[i]

  converted_values = model.convertCubeParameters(
    cube_parameters)
  
  global parameters
  parameters = converted_values[0]

  global physical_parameters
  physical_parameters = converted_values[1]

  for i in range(nparams):
    cube[i] = parameters[i]
  
  return None



def loglike(cube, ndim, nparams):
  
  log_like = model.computeLikelihood(physical_parameters)

  return log_like


print("Starting retrieval\n")
pymultinest.run(
  loglike, 
  priors, 
  model.nbParameters(), 
	resume = False, 
  verbose = model_config.multinest_feedback, 
  importance_nested_sampling = True, 
  sampling_efficiency = model_config.multinest_efficiency, 
  n_live_points = model_config.multinest_nb_living_points, 
  max_iter = model_config.multinest_nb_iterations,
  outputfiles_basename=multinest_output_folder)


#Define the post process configuration
save_post_temperatures = False
save_post_spectra = True
save_post_chemistry = []

postprocess_config = pybear.TransmissionPostProcessConfig(
  save_post_temperatures,
  save_post_spectra,
  save_post_chemistry)


#create a pyBeAR retrieval post process object
post_process = pybear.PostProcess(
  model_config, 
  forward_model_config, 
  postprocess_config,
  observations,
  priors_config)

print("Starting post process\n")
post_process.run(multinest_output_folder + "post_equal_weights.dat")
