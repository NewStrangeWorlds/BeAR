import os
import sys

current_directory = os.path.dirname(os.path.realpath(__file__))
parent_directory = os.path.dirname(current_directory)
sys.path.append(parent_directory)

from lib import pybear
import numpy as np
import pymultinest

#For simplicity, we use the configuration files from the example
#We could also set up the configuration manually as shown in 
#the pyBear_retrieval_3.py example

#setting the basic properties of the model
retrieval_folder = "FlatLineExample/"

#load the retrieval configuration file
model_config = pybear.Config(retrieval_folder)

#create a pyBeAR retrieval post process object
post_process = pybear.PostProcess(model_config)


#we manually load the posteriors from the example
posterior_data = np.loadtxt("FlatLineExample/post_equal_weights.dat")
posteriors = posterior_data[:, 0:-1]
log_like = posterior_data[:, -1]

nb_samples = posteriors.shape[0]

#we can obtain all spectra, including the high-resolution spectra
return_high_res_spectrum = True

#the corresponding wavelength grid:
wavelengths_high_res = np.array(post_process.spectral_grid.wavelength_list)

for i in range(nb_samples):
  physical_parameters = post_process.convertToPhysicalParameters(posteriors[i])

  output = post_process.computeModel(physical_parameters, return_high_res_spectrum)
  
  #the spectra binned to the corresponding observations would be:
  for j in range(len(output.spectrum_obs)):
    model_spectrum_obs = np.array(output.spectrum_obs[j])

  #the high-resolution spectrum is
  spectrum_high_res = np.array(output.spectrum)
  