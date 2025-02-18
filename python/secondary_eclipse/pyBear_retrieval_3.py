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
model_type = "secondary_eclipse"
cross_section_file_path = "/media/data/opacity_data/helios-k/"
spectral_discretisation = 'const_wavenumber'
resolution = 1.0

retrieval_folder = "SecondaryEclipseExample/"
multinest_output_folder = "SecondaryEclipseExample/"
post_output_folder = "SecondaryEclipseExample/"

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
model_config.multinest_efficiency = 0.8 
model_config.multinest_nb_living_points = 800
model_config.multinest_nb_iterations = 0
model_config.multinest_feedback = True
model_config.nb_omp_processes = nb_omp_threads


#read in the observation data and the filter response functions where requried
obs1 = np.loadtxt(retrieval_folder+"wasp-43b_wfc3.dat", skiprows=11)

obs2 = np.atleast_2d(np.loadtxt(retrieval_folder+"wasp-43b_grond_k.dat", skiprows=11))
resp2 = np.loadtxt(retrieval_folder+"../telescope_data/LaSilla_GROND_K_bandpass.dat", skiprows=6)

obs3 = np.atleast_2d(np.loadtxt(retrieval_folder+"wasp-43b_hawki_nb1190.dat", skiprows=11))
resp3 = np.loadtxt(retrieval_folder+"../telescope_data/Paranal_HAWKI_NB1190_bandpass.dat", skiprows=6)

obs4 = np.atleast_2d(np.loadtxt(retrieval_folder+"wasp-43b_hawki_nb2090.dat", skiprows=11))
resp4 = np.loadtxt(retrieval_folder+"../telescope_data/Paranal_HAWKI_NB2090_bandpass.dat", skiprows=6)

obs5 = np.atleast_2d(np.loadtxt(retrieval_folder+"wasp-43b_wircam_h.dat", skiprows=11))
resp5 = np.loadtxt(retrieval_folder+"../telescope_data/CFHT_Wircam_H_bandpass.dat", skiprows=6)

obs6 = np.atleast_2d(np.loadtxt(retrieval_folder+"wasp-43b_wircam_ks.dat", skiprows=11))
resp6 = np.loadtxt(retrieval_folder+"../telescope_data/CFHT_Wircam_Ks_bandpass.dat", skiprows=6)

obs7 = np.atleast_2d(np.loadtxt(retrieval_folder+"wasp-43b_spitzer_1.dat", skiprows=11))
resp7 = np.loadtxt(retrieval_folder+"../telescope_data/Spitzer_irac1_bandpass.dat", skiprows=6)

obs8 = np.atleast_2d(np.loadtxt(retrieval_folder+"wasp-43b_spitzer_2.dat", skiprows=11))
resp8 = np.loadtxt(retrieval_folder+"../telescope_data/Spitzer_irac2_bandpass.dat", skiprows=6)

#and create the pyBeAR observation inputs
wfc3_obs = pybear.Observation(
  "WFC3", "band-spectroscopy", obs1[:,0:2], obs1[:,2], obs1[:,3])

grond_k_obs = pybear.Observation(
  "GROND_K", "photometry", obs2[:,0:2], obs2[:,2], obs2[:,3])
grond_k_obs.filter_response = np.transpose(resp2)
grond_k_obs.filter_detector_type = "energy"

hawki_1_obs = pybear.Observation(
  "HAWKI_NB1190", "photometry", obs3[:,0:2], obs3[:,2], obs3[:,3])
hawki_1_obs.filter_response = np.transpose(resp3)
hawki_1_obs.filter_detector_type = "energy"

hawki_2_obs = pybear.Observation(
  "HAWKI_NB2090", "photometry", obs4[:,0:2], obs4[:,2], obs4[:,3])
hawki_2_obs.filter_response = np.transpose(resp4)
hawki_2_obs.filter_detector_type = "energy"

wircam_h_obs = pybear.Observation(
  "Wircam_H", "photometry", obs5[:,0:2], obs5[:,2], obs5[:,3])
wircam_h_obs.filter_response = np.transpose(resp5)
wircam_h_obs.filter_detector_type = "energy"

wircam_ks_obs = pybear.Observation(
  "Wircam_Ks", "photometry", obs6[:,0:2], obs6[:,2], obs6[:,3])
wircam_ks_obs.filter_response = np.transpose(resp6)
wircam_ks_obs.filter_detector_type = "energy"

spitzer_1_obs = pybear.Observation(
  "Spitzer_1", "photometry", obs7[:,0:2], obs7[:,2], obs7[:,3])
spitzer_1_obs.filter_response = np.transpose(resp7)
spitzer_1_obs.filter_detector_type = "photon"

spitzer_2_obs = pybear.Observation(
  "Spitzer_2", "photometry", obs8[:,0:2], obs8[:,2], obs8[:,3])
spitzer_2_obs.filter_response = np.transpose(resp8)
spitzer_2_obs.filter_detector_type = "photon"


#create a list of all observations
observations = list([
  wfc3_obs,
  grond_k_obs,
  hawki_1_obs,
  hawki_2_obs,
  wircam_h_obs,
  wircam_ks_obs,
  spitzer_1_obs,
  spitzer_2_obs])


#create the list of priors
priors_config = list([
  pybear.Prior("delta", "log_g", [3.64464]),
  pybear.Prior("delta", "Rp/Rs", [0.151271]),
  pybear.Prior("uniform", "M/H", [0.1, 2.45]),
  pybear.Prior("delta", "C/O", [0.55]),
  pybear.Prior("uniform", "temperature1", [5000, 1000]),
  pybear.Prior("uniform", "temperature2", [0.1, 1.0]),
  pybear.Prior("uniform", "temperature3", [0.1, 2.0]),
  pybear.Prior("uniform", "temperature4", [0.1, 2.0]),
  pybear.Prior("uniform", "temperature5", [0.1, 2.0]),
  pybear.Prior("uniform", "temperature6", [0.1, 2.0]),
  pybear.Prior("uniform", "temperature7", [0.1, 2.0])])


#Now we configure the forward model
nb_grid_points = 70
bottom_pressure = 100
top_pressure = 1e-3

opacity_species_data = np.array([
  ['CIA-H2-H2', 'CIA/H2-H2'],
  ['CIA-H2-He', 'CIA/H2-He'],
  ['Na',        'Alkali_Allard/Na'],
  ['K',         'Alkali_Allard/K'],
  ['H2O',       'Molecules/H2O_pokazatel'],
  ['CO',        'Molecules/12C-16O__Li2015_e2b'],
  ['TiO',       'Molecules/48Ti-16O__Toto_e2b_v2'],
  ['VO',        'Molecules/51V-16O__VOMYT_e2b'],
  ['SH',        'Molecules/32S-1H__GYT_e2b'],
  ['H2S',       'Molecules/1H2-32S__AYT2_e2b'],
  ['FeH',	      'Molecules/56Fe-1H__Yueqi_e2b'],
  ['CH4',       'Molecules/12C-1H4__YT34to10_e2b'],
  ['CO2',       'Molecules/12C-16O2__CDSD_4000_e2b'],
  ['HCN',       'Molecules/1H-12C-14N__Harris_e2b'],
  ['MgH',       'Molecules/26Mg-1H__Yadin_e2b'],
  ['TiH',       'Molecules/48Ti-1H__MoLLIST_e2b'],
  ['CrH',       'Molecules/52Cr-1H__Yueqi_e2b'],
  ['CaH',       'Molecules/40Ca-1H__MoLLIST_e2b']])
  
#configuration of the temperature profile
temperature_profile = "poly"
#a constant temperature profile does not need any config parameters
temperature_profile_parameters = ["6", "1"] 

#radiative transfer model
rt_model = "scm"
rt_model_parameters = [""]

#configuration of the chemistry
#here, we use the equilibrium chemistry
chemistry_models = ["eq"]

#the chemical species for each of the two models
chemistry_parameters = [[retrieval_folder + "fastchem_parameters.dat"]]

#stellar spectrum model
stellar_model = "file"
stellar_model_parameters = ["SecondaryEclipseExample/WASP-43.dat"]

#cloud model
cloud_models = []
cloud_models_parameters = [[""]]

opacity_species_symbols = opacity_species_data[:, 0]
opacity_species_folders = opacity_species_data[:, 1]

#Create the configuration of the forward model
forward_model_config = pybear.SecondaryEclipseConfig(
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
  opacity_species_folders,
  stellar_model,
  stellar_model_parameters,
  cloud_models,
  cloud_models_parameters)


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
save_post_spectra = True
save_contribution_functions = True
save_post_chemistry = ['H2O', 'CO', 'CO2', 'CH4', 'NH3']

postprocess_config = pybear.SecondaryEclipsePostProcessConfig(
  save_post_temperatures,
  save_post_spectra,
  save_contribution_functions,
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
