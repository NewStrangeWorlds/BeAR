Setting up a retrieval
======================

To start a retrieval calculation, BeAR requires the following files in the folder the executable is called with:

  - ``retrieval.config`` - the main configuration file for the retrieval
  
  - ``forward_model.config`` - the configuration file of chosen forward forward_model
  
  - ``priors.config`` - the setup list for the prior distributions of the free parameters
  
  - ``observations.list`` - the list of observational data files that the retrieval should use

  
Main retrieval file
...................

The file ``retrieval.config`` contains the basic information for the retrieval setup. 

.. include:: ../examples/retrieval.config_example
   :literal:
   
The following parameters need to be set:

| ``Use GPU``
|  Determines the use of the graphics card. Set either ``Y`` or ``1`` to run the calculation on the GPU. 
  Any other input is interpreted as running purely on the CPU.
    
| ``OpenMP processor number`` 
|  Sets the number of processor cores used for parallel computing on the CPU. Some parts
  of BeAR will still run on the CPU, even if ``Use GPU`` is enabled. This, for example, is the case for the
  FastChem chemistry code that doesn't run GPUs. BeAR will run certain calculations in parallel on the CPU as well, using
  OpenMP. Set this parameter to ``0`` if you want to use all available cores. Note that OpenMP can only use a single, 
  multi-core processor.
    
| ``forward model type`` 
|  Sets the forward model that is supppsed to be used. BeAR currently supports the following models:

     - ``transmission`` - Transmission spectrum
    
     - ``secondary_eclipse`` - Secondary eclipse / occulation spectrum
     
     - ``emission`` - Emission spectrum
     
     - ``flat_line`` - Fits a flat line to the data
  
  Descriptions of the forward models can be found :ref:`here <sec:forward_models>` 

| ``Spectral grid parametrisation``
|  Sets the parametrisation of the spectral grid that will be used for the computation of the 
  high-resolution spectrum. This high-resolution grid should generally be finer than that of the
  observational data.. The following options are available:

     - ``const_wavelength x`` - a constant step in wavelength space with a step size of ``x`` in :math:`\mathrm{\mu m}`
     
     - ``const_wavenumber x`` - a constant step in wavenumber space with a step size of ``x`` :math:`\mathrm{cm}^{-1}`
     
     - ``const_resolution x`` - a constant spectral resolution :math:`x = \lambda/\Delta\lambda`

| ``Opacity data folder``
|  Location of the folder with opacities for the gas species. Details on the required format
    of the opacity data can be found in this :ref:`section <sec:opacity_data>` 
    
| ``Use error inflation prior``
|  Determines the use of the error inflation. This will artifically enlarge the error bars of the
  observational data and, thus, will generally make it easier for the retrieval to find a
  solution. The use of the error inflation acknowledges that certain physical or chemical processes
  are missing from the simple forward model of the retrieval. The form of the employed error inflation
  is described in `Kitzmann et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020ApJ...890..174K/>`_

The remaining parameters refer to the nested-sampling code MultiNest. For a description of the MultiNest code
and its parameters, we refer to `Feroz & Hobs (2009) <https://ui.adsabs.harvard.edu/abs/2008MNRAS.384..449F/>`_
and `Feroz et al. (2008) <https://ui.adsabs.harvard.edu/abs/2009MNRAS.398.1601F/>`_ .

The parameters that can be set here include:

| ``Importance nested sampling``
|  Turns the use of importance nested sampling on or off. To use importance nested sampling, set this parameter 
  to either ``Y`` or ``1``. Any other value is interpreted as ``N``. Importance nested sampling requires
  a bit more memory but, on the other hand, also increases the convergence speed and
  the overall accuracy of the Bayesian evidence calculation. Unless memory is a real
  bottleneck, there should be no reason to deactivate important nested sampling.
  
| ``Mode separation``
|  MultiNest has the ability to trace different modes in a posterior distribution and to save them separately. 
  This option turns the use of mode separation on or off. To use it, set this parameter 
  to either ``Y`` or ``1``. Any other value is interpreted as ``N``. Note that this option in untested and
  might not work properly.
  
| ``Number of live points``
|  Sets the number of live points uses by MultiNest. Generally speaking, a high-dimensional parameter space
  requires a higher amount of live points. It is strongly recommended to perform sensitivy tests by
  increasing this number and check if the posterior distributions have converged.
  
| ``Efficiency``
|  Sets the efficiency that determines the way MultiNest draws new points from the parameter space. 
  For more details on this parameter check the MultiNest documentation. The authors of MultiNest suggest 
  to use an efficiency of ``0.8`` for parameter estimations and ``0.3`` when the Bayesian evidence is 
  wanted at a high accuracy.
  
| ``Maximum number of iterations``
|  The maximum number of iterations MultiNest will use before the nested sampling is stopped. 
  A value of ``0`` indicates that MultiNest will perform the nested sampling until its convergence criteria are met.
  
| ``Resume``
|  If this is set to ``Y`` or ``1``, MultiNest will try to resume a previously started retrieval run. The files MultiNest 
  needs to restart the nested sampling must all be present in the retrieval folder. This option is usefull if BeAR 
  is run on a cluster with a strict time limit. A restart can also be used if a previous MultiNest run was stopped 
  at its maximum number of iterations.
  
| ``Console feedback``
|  MultiNest will regularly report the current total number of model evaluations and estimates for the Bayesian evidences 
  when this option is turned on with ``Y`` or ``1``.
  
| ``Print parameter values and likelihoods``
|  Determines whether the parameter values and the computed likelihood values for all models should be displayed 
  (``Y`` or ``1``).  If BeAR is run on a cluster and the terminal output is redirected to a file, 
  it is usually a good idea to deactivate this option. Otherwise, the output file could become quite large.


Forward model configuration file
................................

The file ``forward_model.config`` contains the configuration for the forward model.
Its structure depends on the chosen model and is discussed in the :ref:`section <sec:forward_models>` on forward models.


Prior distributions file
........................

The ``priors.config`` file contains the information on the prior distributions of the free parameters.
The file has the following structure:

.. include:: ../examples/priors_example1.config
   :literal:

The first column lists the distribution type of the prior, the second column the model parameter name, and the 
remaining columns the parameters of the distribution.

BeAR supports the following types of distributions:

  - ``delta`` - a delta distribution with a fixed value

  - ``uniform`` - a uniform distribution with a lower and upper bound

  - ``log_uniform`` - a log-uniform distribution with a lower and upper bound

  - ``gaussian`` - a Gaussian distribution with a mean and standard deviation

  - ``linked`` - links this prior to that of another parameter

The type, order and number of these parameters in prior distributions file is determined by chosen forward model.

A special case is the ``linked`` distribution. This distribution links the prior distribution 
of one parameter to that of another. The config parameter for this distribution is the line number of the
parameter distribution that it should be linked to.
An example of this is shown below:

.. include:: ../examples/priors_example2.config
   :literal:

Here, the prior for the CO2 mixing ratios is linked to the fifth parameter, which is the mixing ratio of CH4 in this example.
Thus, CO2 will always have the same mixing ratio as methane for this retrieval setup. It is important to note that BeAR cannot 
check the consistency of the linked parameters. For example, if the linked parameter is a temperature, the resulting mixing ratios
of CO2 would make no sense. It is the user's responsibility to ensure that the linked parameters are consistent.

Currently, BeAR assumes the following units for the free model parameters:
  
  - ``log g`` - logarithm of the surface gravity in :math:`\mathrm{cm/s^2}`

  - ``R_p`` - planetary radius in :math:`\mathrm{R_{Jup}}`

  - ``R_s`` - stellar radius in :math:`\mathrm{R_{Sun}}`

  - ``distance`` - distance of the object in :math:`\mathrm{pc}`


Observational data file
.......................

The ``observations.list`` file contains a list of data files with the observational data that the retrieval should use.
An example is shown below:

.. include:: ../examples/observations_example.list
   :literal:

BeAR can use multiple observational data files at the same time. The observations do not need to be orderer in any specific way.
They also do not need to be continuous in wavelength space, gaps are are allowed to be present between the different observations.
It also possible to mix different observational types, e.g. photometric data together with spectroscopic data. The format of the 
these files is described in the :ref:`section <sec:observational_data>` on observational data.