import os
import sys

current_directory = os.path.dirname(os.path.realpath(__file__))
parent_directory = os.path.dirname(current_directory)
sys.path.append(parent_directory)

from lib import bear


#setting the basic properties of the model
retrieval_folder = "SecondaryEclipseExampleBB/"

#load the retrieval configuration file
model_config = bear.Config(retrieval_folder)

#create a pyBeAR retrieval object
model = bear.Retrieval(model_config)

print("Starting retrieval\n")
model.run()

#create a pyBeAR retrieval post process object
post_process = bear.PostProcess(model_config)

print("Starting post process\n")
post_process.run()
