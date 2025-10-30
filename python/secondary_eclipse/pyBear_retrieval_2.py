import os
import sys

current_directory = os.path.dirname(os.path.realpath(__file__))
parent_directory = os.path.dirname(current_directory)
sys.path.append(parent_directory)

from lib import pybear
import numpy as np
import pymultinest


#setting the basic properties of the model
retrieval_folder = "SecondaryEclipseExample/"

#load the retrieval configuration file
model_config = pybear.Config(retrieval_folder)

#create a pyBeAR retrieval object
model = pybear.Retrieval(model_config)


#load the observation data to compute the likelihood
obs1 = np.loadtxt(retrieval_folder+"wasp-43b_wfc3.dat", skiprows=11)
obs2 = np.atleast_2d(np.loadtxt(retrieval_folder+"wasp-43b_grond_k.dat", skiprows=11))
obs3 = np.atleast_2d(np.loadtxt(retrieval_folder+"wasp-43b_hawki_nb1190.dat", skiprows=11))
obs4 = np.atleast_2d(np.loadtxt(retrieval_folder+"wasp-43b_hawki_nb2090.dat", skiprows=11))
obs5 = np.atleast_2d(np.loadtxt(retrieval_folder+"wasp-43b_wircam_h.dat", skiprows=11))
obs6 = np.atleast_2d(np.loadtxt(retrieval_folder+"wasp-43b_wircam_ks.dat", skiprows=11))
obs7 = np.atleast_2d(np.loadtxt(retrieval_folder+"wasp-43b_spitzer_1.dat", skiprows=11))
obs8 = np.atleast_2d(np.loadtxt(retrieval_folder+"wasp-43b_spitzer_2.dat", skiprows=11))

#create a list of observations and their errors
#we have to use the same order as in the observations.list file
obs_data = list([
  obs1[:,2],
  obs2[:,2],
  obs3[:,2],
  obs4[:,2],
  obs5[:,2],
  obs6[:,2],
  obs7[:,2],
  obs8[:,2]])

obs_error = list([
  obs1[:,3],
  obs2[:,3],
  obs3[:,3],
  obs4[:,3],
  obs5[:,3],
  obs6[:,3],
  obs7[:,3],
  obs8[:,3]])

nb_obs = len(obs_data)

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
  
  log_like = 0
  
  for i in range(nb_obs):
    model_spectrum_obs = np.array(output.spectrum_obs[i])
    error = obs_error[i]
    obs_delta = obs_data[i] - model_spectrum_obs

    for j in range(obs_delta.size):
      error_square = error[j]**2
      log_like += (- 0.5 * np.log(error_square* 2.0 * np.pi)
                   - 0.5 * obs_delta[j]**2 / error_square)

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
