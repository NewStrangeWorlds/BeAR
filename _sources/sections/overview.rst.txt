Overview
~~~~~~~~

``BeAR`` is an open-source computer program that can perform retrieval calculations for observational 
data of planetary atmospheres as well as brown dwarfs. It is the successor of the previous retrieval 
code ``HELIOS-r2`` (`Kitzmann et al. 2020 <http://adsabs.harvard.edu/abs/2020ApJ...890..174K>`_) that 
has been successfully applied for the characterisation of data from ground-based and space telescopes 
for a wide range of objects.

``BeAR`` is currently able to perform retrievals for emission, transmission, and secondary-eclipse spectra. 
It can use low to mid-resolution spectral data, as well as photometric data points. It has support for 
constant temperature and abundance profiles of atmospheric constituents but can also retrieve 
vertically-varying profiles if the observational data is good enough.

While the core of ``BeAR`` is written in C++ and CUDA, it also comes with a Python interface, ``pyBeAR``, that 
allows the forward model and retrieval calculations to be performed in Python. The interface allows ``BeAR`` to be
coupled to parameter space explorers or machine-learning algorithms written in Python. The stand-alone C++ version
is currently only able to use the MultiNest sampler for the retrieval calculations, while the Python interface
also supports the use of other sampling algorithms, such as Dynesty or MCMC.
