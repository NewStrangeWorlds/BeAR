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
model_type = "flat_line"
cross_section_file_path = ""
spectral_discretisation = 'const_wavenumber'
resolution = 0.5

retrieval_folder = "FlatLineExample/"
multinest_output_folder = "FlatLineExample/"
post_output_folder = "FlatLineExample/"

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
obs = np.loadtxt(retrieval_folder+"55Cnce_V5.dat", skiprows=11)

nircam = pybear.Observation(
  "NIRCAM", "spectroscopy", obs[:,0], obs[:,1], obs[:,2])

#create a list of all observations
observations = list([nircam])


#create the list of priors
priors_config = list([
  pybear.Prior("uniform", "eclipse_depth", [0, 200])])


#Create the configuration of the forward model
#The flat line model does not require any additional configuration
forward_model_config = None


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
save_post_spectra = True

postprocess_config = pybear.FlatLinePostProcessConfig(
  save_post_spectra)


#create a pyBeAR retrieval post process object
post_process = pybear.PostProcess(
  model_config, 
  forward_model_config, 
  postprocess_config,
  observations,
  priors_config)

print("Starting post process\n")
post_process.run(multinest_output_folder + "post_equal_weights.dat")
