# BeAR - The Bern Atmospheric Retrieval code
#### Authors: Daniel Kitzmann ####


# Overview #
BeAR is an open source model that can perform atmospheric retrieval of brown dwarf and exoplanet spectra. It is based on the previous retrieval code HELIOS-r2 that has been introduced and described in Kitzmann et al. (2020).

BeAR uses a Bayesian statistics approach by employing a nested sampling method to generate posterior distributions and calculate the Bayesian evidence. The nested sampling itself is done by the Multinest library (https://github.com/farhanferoz/MultiNest). The computationally most demanding parts of the model have been written in NVIDIA's CUDA language for an increase in computational speed. BeAR can work on both, pure CPU as well as hybrid CPU/GPU setups. Running it purely on a CPU is not recommended, though, as the runtimes can be  by a factor of 10 or 100 longer compared to running it on a GPU.

Successful applications include retrieval of brown dwarf emission spectra (Kitzmann et al. 2020) and secondary eclipse measurements of exoplanets (Bourrier et al. 2020).


# User Guide #

BeAR comes with a user guide that can be found here: https://newstrangeworlds.github.io/BeAR/ . It describes the installation and usage of BeAR, including a details on the different forward models and additional modules it contains. The documentation is currently under heavy revision and might not yet be fully complete.


# Python Interface #

Since its newest version, BeAR now includes a Python interface called pyBeAR that allows retrieval calculations to be done from within Python. Examples for the different forward models are provided in the python folder. In addition to performing retrievals, the Python interface can also call the forward models directly. That way, for example, transmission or emission spectra can be calculated for user-defined temperature profiles or chemical compositions. A detailed description of pyBeAR is provided in the user guide.
