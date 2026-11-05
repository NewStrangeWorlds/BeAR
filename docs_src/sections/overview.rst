Overview
~~~~~~~~

``BeAR`` is an open-source computer program that can perform retrieval calculations for observational 
data of planetary atmospheres as well as brown dwarfs. It is the successor of the previous retrieval 
code ``HELIOS-r2`` (`Kitzmann et al. 2020 <http://adsabs.harvard.edu/abs/2020ApJ...890..174K>`_) that 
has been successfully applied for the characterisation of data from ground-based and space telescopes 
for a wide range of objects.

``BeAR`` is currently able to perform retrievals for emission, transmission, secondary-eclipse, and
high-resolution phase-curve spectra.
It can use low to mid-resolution spectral data, as well as photometric data points. It has support for
constant temperature and abundance profiles of atmospheric constituents but can also retrieve
vertically-varying profiles if the observational data is good enough.

The core of ``BeAR`` is written in C++ and CUDA, but the code is used through its Python module
``bear``, which allows the forward model and retrieval calculations to be performed in Python. The
interface allows ``BeAR`` to be coupled to parameter space explorers or machine-learning algorithms
written in Python. Retrievals are driven from Python by a sampler such as MultiNest (via PyMultiNest),
nautilus, or dynesty. MultiNest can optionally be run with MPI to distribute the likelihood evaluations
over several GPUs.
