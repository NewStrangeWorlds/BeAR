import os
import sys

current_directory = os.path.dirname(os.path.realpath(__file__))
parent_directory = os.path.dirname(current_directory)
sys.path.append(parent_directory)

from lib import pybear
import numpy as np
import pymultinest


#setting the basic properties of the model
retrieval_folder = "../../TransmissionExample/"

#load the retrieval configuration file
model_config = pybear.Config(retrieval_folder)

#create a pyBeAR retrieval object
model = pybear.Retrieval(model_config)


#load the observation data to compute the likelihood
observation = np.loadtxt(retrieval_folder+"WASP-12b_kreidberg.dat", skiprows=11)
obs_data = observation[:,2]
obs_error = observation[:,3]

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


#Instead of using the computeLikelihood function from BeAR, 
#we calculate the likelihood using spectra provided by BeAR and 
#the corresponding observations
def loglike(cube, ndim, nparams):

  output = model.computeModel(physical_parameters, False)
  
  model_spectrum_obs = np.array(output.spectrum_obs[0])

  obs_delta = obs_data - model_spectrum_obs

  log_like = 0

  for i in range(obs_data.size):
    error_square = obs_error[i]**2
    log_like += (- 0.5 * np.log(error_square* 2.0 * np.pi)
                 - 0.5 * obs_delta[i]**2 / error_square)

  return log_like



print("Starting retrieval\n")
pymultinest.run(
  loglike, 
  priors, 
  model.nbParameters(), 
	resume = False, 
  verbose = True, 
  importance_nested_sampling = True, 
  sampling_efficiency = model_config.multinest_efficiency, 
  n_live_points = model_config.multinest_nb_living_points, 
  max_iter = model_config.multinest_nb_iterations,
  outputfiles_basename=retrieval_folder)


#create a pyBeAR retrieval post process object
post_process = pybear.PostProcess(model_config)

print("Starting post process\n")
post_process.run()
