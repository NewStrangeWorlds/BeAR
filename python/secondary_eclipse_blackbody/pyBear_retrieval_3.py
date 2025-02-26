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
model_type = "secondary_eclipse_bb"
cross_section_file_path = ""
spectral_discretisation = 'const_wavenumber'
resolution = 1.0

retrieval_folder = "SecondaryEclipseExampleBB/"
multinest_output_folder = "SecondaryEclipseExampleBB/"
post_output_folder = "SecondaryEclipseExampleBB/"

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
obs = np.atleast_2d(np.loadtxt(retrieval_folder+"MIRI.dat", skiprows=11))

miri = pybear.Observation(
  "MIRI", "photometry", obs[:,0:2], obs[:,2], obs[:,3])
miri_resp = np.loadtxt(retrieval_folder+"../telescope_data/JWST_MIRI_F1500W.dat", skiprows=6)
miri.filter_response = np.transpose(miri_resp)
miri.filter_detector_type = "photon"


#create a list of all observations
observations = list([miri])


#create the list of priors
priors_config = list([
  pybear.Prior("uniform", "planet_temperature", [500, 3000]),
  pybear.Prior("gaussian", "radius_ratio", [0.03385, 0.001488]),
  pybear.Prior("gaussian", "stellar_temperature", [3496.0, 25.0])])

#stellar spectrum model
stellar_model = "blackbody"
stellar_model_parameters = [] #blackbody does not require any parameters


#Create the configuration of the forward model
forward_model_config = pybear.OccultationBlackBodyConfig(
  stellar_model,
  stellar_model_parameters)


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

postprocess_config = pybear.OccultationBlackBodyPostConfig(
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
