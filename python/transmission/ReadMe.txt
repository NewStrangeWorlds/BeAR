This folder contains several Python scripts that show how use the *transmission spectroscopy* model and retrieval with pyBeAR.

The following examples are provided:

- pyBear_transmission_spectrum.py
This Python script shows how to call the transmission spectroscopy forward model directly for a given p-T structure, chemical composition, and cloud layer.

- pyBear_retrieval_0.py
This script runs the TransmissionExample retrieval case in a stand-alone fashion. It is the Python equivalent of the C++ version of BeAR.

- pyBear_retrieval_1.py
This script runs the TransmissionExample retrieval case by using the pyMultiNest module, the Python version of MultiNest. It shows how the prior and likelihood functions of BeAR are used and interfaced with pyMultiNest.

- pyBear_retrieval_2.py
This script runs the TransmissionExample retrieval case by using again pyMultiNest. In addition to the previous script, the likelihood is here directly calculated within Python.

- pyBear_retrieval_3.py
This script shows how to use pyBeAR without any of the normal configuration files. The entire retrieval and model configuration is done in Python. The parameters are equivalent to those found in the configuration files of the TransmissionExample folder.

- pyBear_postprocess.py
This example shows how to obtain spectra and atmospheric structures to perform a custom post-process.
