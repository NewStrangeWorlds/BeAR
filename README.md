# Helios-r2 - A Bayesian Nested-Sampling Retrieval Code
#### Authors: Daniel Kitzmann ####

# Overview #
Helios-r2 is an open source model that can perform atmospheric retrieval of brown dwarf and exoplanet spectra. It has been introduced and described in Kitzmann et al. (2020).
and is the successor to the original Helios-r code described by Lavie et al (2017). This original version, however, has never been publicly released. The new version has been completely written from scratch in C++/CUDA and includes various improvements over the original one.

Helios-r2 uses a Bayesian statistics approach by employing a nested sampling method to generate posterior distributions and calculate the Bayesian evidence. The nested sampling itself is done by the Multinest library (https://github.com/farhanferoz/MultiNest). The computationally most demanding parts of the model have been written in NVIDIA's CUDA language for an increase in computational speed. Helios-r2 can work on both, pure CPU as well as hybrid CPU/GPU setups. Running it purely on a CPU is not recommended, though, as the runtimes can be  by a factor of 10 or 100 longer compared to running it on a GPU.

Successful applications include retrieval of brown dwarf emission spectra (Kitzmann et al. 2020) and secondary eclipse measurements of exoplanets (Bourrier et al. 2020).

A complete guide to the model can be found in the file **helios-r2_guide.pdf** located within this repository.
