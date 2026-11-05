This folder contains several Python scripts that show how use the *secondary eclipse blackbody* model and retrieval with BeAR. 
The relative paths in these files are defined for the Python scripts to be run from the BeAR root directory rather than from the local Python from.

The following examples are provided:

- bear_retrieval_0.py
This script runs the SecondaryEclipseExampleBB retrieval case in a stand-alone fashion. It is the Python equivalent of the C++ version of BeAR.

- bear_retrieval_1.py
This script runs the SecondaryEclipseExampleBB retrieval case by using the pyMultiNest module, the Python version of MultiNest. It shows how the prior and likelihood functions of BeAR are used and interfaced with pyMultiNest.

- bear_retrieval_2.py
This script runs the SecondaryEclipseExampleBB retrieval case by using again pyMultiNest. In addition to the previous script, the likelihood is here directly calculated within Python.

- bear_retrieval_3.py
This script shows how to use BeAR without any of the normal configuration files. The entire retrieval and model configuration is done in Python. The parameters are equivalent to those found in the configuration files of the SecondaryEclipseExampleBB folder.

- bear_postprocess.py
This example shows how to obtain posterior spectra to perform a custom post-process.
