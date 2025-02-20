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
nb_omp_threads = 0
model_type = "emission"
cross_section_file_path = "/media/data/opacity_data/helios-k/"
spectral_discretisation = 'const_wavenumber'
resolution = 1.0

retrieval_folder = "EmissionExample/"
multinest_output_folder = "EmissionExample/"
post_output_folder = "EmissionExample/"


#create the general model config
model_config = pybear.Config(
  use_gpu, 
  model_type, 
  cross_section_file_path, 
  spectral_discretisation, 
  resolution,
  multinest_output_folder,
  post_output_folder)

#configure additional parameters
model_config.use_error_inflation = True
model_config.multinest_efficiency = 0.8 
model_config.multinest_nb_living_points = 4000
model_config.multinest_nb_iterations = 0
model_config.multinest_feedback = True
model_config.nb_omp_processes = nb_omp_threads


#read in the observation data
obs = np.loadtxt(retrieval_folder+"gj570d_spex.dat", skiprows=11)
obs_wavelengths = obs[:,0]
obs_data = obs[:,1]
obs_error = obs[:,2]
obs_line_spread = obs[:,3]
log_like_weight = obs[:,4]

#and create the pyBeAR observation inputs
spex = pybear.Observation(
  "GJ570D SpeX", "spectroscopy", obs_wavelengths, obs_data, obs_error)

spex.likelihood_weight = log_like_weight
spex.instrument_profile_fwhm = obs_line_spread

#create a list of all observations
observations = list([spex])


#create the list of priors
priors_config = list([
  pybear.Prior("uniform", "log_g", [4.0, 6.0]),
  pybear.Prior("uniform", "scaling_factor", [0.1, 5.0]),
  pybear.Prior("gaussian", "distance", [5.8819, 0.0029], "pc"),
  pybear.Prior("log_uniform", "MR_H2O", [1e-12, 0.01]),
  pybear.Prior("log_uniform", "MR_CH4", [1e-12, 0.01]),
  pybear.Prior("log_uniform", "MR_NH3", [1e-12, 0.01]),
  pybear.Prior("log_uniform", "MR_K", [1e-12, 0.01]),
  pybear.Prior("log_uniform", "MR_H2S", [1e-12, 0.01]),
  pybear.Prior("log_uniform", "MR_CO2", [1e-12, 0.01]),
  pybear.Prior("log_uniform", "MR_CO", [1e-12, 0.01]),
  pybear.Prior("uniform", "temperature1", [5000, 1000]),
  pybear.Prior("uniform", "temperature2", [0.3, 0.95]),
  pybear.Prior("uniform", "temperature3", [0.3, 0.95]),
  pybear.Prior("uniform", "temperature4", [0.4, 0.95]),
  pybear.Prior("uniform", "temperature5", [0.4, 0.95]),
  pybear.Prior("uniform", "temperature6", [0.4, 0.95]),
  pybear.Prior("uniform", "temperature7", [0.4, 0.95])])

#since we use error inflation here, we need to add an additional prior
max_error = np.max(obs_error)
min_error = np.min(obs_error)

max_error = np.log10(100.0 * max_error * max_error)
min_error = np.log10(0.1 * min_error * min_error)

priors_config.append(
  pybear.Prior("uniform", "error_exponent", [min_error, max_error]))


#Now we configure the forward model
nb_grid_points = 70
bottom_pressure = 300
top_pressure = 1e-3

opacity_species_data = np.array([
  ['CIA-H2-H2', 'CIA/H2-H2'],
  ['CIA-H2-He', 'CIA/H2-He'],
  ['Na',        'Alkali_Allard/Na'],
  ['K',         'Alkali_Allard/K'],
  ['H2O',       'Molecules/1H2-16O__POKAZATEL_e2b_n'],
  ['CH4',       'Molecules/12C-1H4__YT34to10_e2b'],
  ['CO2',       'Molecules/12C-16O2__CDSD_4000_e2b'],
  ['CO',        'Molecules/12C-16O__Li2015_e2b_n'],
  ['NH3',       'Molecules/14N-1H3__CoYuTe_e2b'],
  ['H2S',       'Molecules/1H2-32S__AYT2_e2b']])
  
#configuration of the temperature profile
temperature_profile = "poly"
temperature_profile_parameters = ["6", "1"] 

#radiative transfer model
rt_model = "scm"
rt_model_parameters = [""]

#configuration of the chemistry
#here, we use the isoprofile and background models
chemistry_models = ["iso", "bg"]

#the chemical species for each of the two models
chemistry_parameters = [
  ["H2O", "CH4", "NH3", "K", "H2S", "CO2", "CO"],
  ["H2He"]]


opacity_species_symbols = opacity_species_data[:, 0]
opacity_species_folders = opacity_species_data[:, 1]

#Create the configuration of the forward model
forward_model_config = pybear.EmissionModelConfig(
  nb_grid_points,
  bottom_pressure,
  top_pressure,
  temperature_profile,
  temperature_profile_parameters,
  rt_model,
  rt_model_parameters,
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
save_post_temperatures = True
save_effective_temperatures = True
save_post_spectra = True
save_contribution_functions = True
save_post_chemistry = []

postprocess_config = pybear.EmissionPostProcessConfig(
  save_post_temperatures,
  save_effective_temperatures,
  save_post_spectra,
  save_contribution_functions,
  save_post_chemistry)


#in order to compute effective temperatures, we need an
#additional "observations" that provides the retrieval with
#a wider wavelength range than the orginal observation
obs_add = pybear.Observation(
  "Postprocess_Spectrum", "band-spectroscopy", [[0.5, 20.0]], [1.0], [1.0])

observations.append(obs_add)

#create a pyBeAR retrieval post process object
post_process = pybear.PostProcess(
  model_config, 
  forward_model_config, 
  postprocess_config,
  observations,
  priors_config)

print("Starting post process\n")
post_process.run(multinest_output_folder + "post_equal_weights.dat")
