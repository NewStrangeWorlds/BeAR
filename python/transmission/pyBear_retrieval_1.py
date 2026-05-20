import os
import sys

current_directory = os.path.dirname(os.path.realpath(__file__))
parent_directory = os.path.dirname(current_directory)
sys.path.append(parent_directory)

from lib import bear
import numpy as np
from lib.bear_multinest_path import MULTINEST_LIB_DIR
import ctypes as _ctypes
_ctypes.CDLL(MULTINEST_LIB_DIR + '/libmultinest.so', mode=_ctypes.RTLD_GLOBAL)
import pymultinest


#setting the basic properties of the model
retrieval_folder = "TransmissionExample/"

#load the retrieval configuration file
model_config = bear.Config(retrieval_folder)

#create a pyBeAR retrieval object
model = bear.Retrieval(model_config)


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
  verbose = True, 
  importance_nested_sampling = True, 
  sampling_efficiency = 0.8, 
  n_live_points = 800, 
  max_iter = 0,
  outputfiles_basename=retrieval_folder)


#create a pyBeAR retrieval post process object
post_process = bear.PostProcess(model_config)

print("Starting post process\n")
post_process.run()
