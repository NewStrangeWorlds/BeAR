import os
import sys

current_directory = os.path.dirname(os.path.realpath(__file__))
parent_directory = os.path.dirname(current_directory)
sys.path.append(parent_directory)

from lib import pybear
import numpy as np
import pickle
from dynesty import DynamicNestedSampler, NestedSampler


#setting the basic properties of the model
retrieval_folder = "trappist-1e_python/"

#load the retrieval configuration file
model_config = pybear.Config(retrieval_folder)

#create a pyBeAR retrieval object
model = pybear.Retrieval(model_config)


#Define the priors and likelihood functions for MultiNest
#We use the internal functions of BeAR to convert the cube parameters
#and to compute the likelihood
nb_param = model.nbParameters()


def ptform(cube):
  cube_parameters = np.zeros(cube.size)

  for i in range(cube.size):
    cube_parameters[i] = cube[i]

  converted_values = model.convertCubeParameters(
    cube_parameters)
  
  parameters = converted_values[0]
  physical_parameters = converted_values[1]
  
  return physical_parameters



def loglike(cube):
  
  log_like = model.computeLikelihood(cube)
  #print(log_like)
  return log_like



print("Starting retrieval\n")
# dsampler = DynamicNestedSampler(loglike, ptform, nb_param, sample='rslice')
# dsampler.run_nested(dlogz_init=0.05, nlive_init=4000, nlive_batch=100, checkpoint_file='dynesty.save')


dsampler = NestedSampler(loglike, ptform, nb_param, sample='rslice', nlive=8000)
dsampler.run_nested(checkpoint_file='dynesty_trappist.save')
