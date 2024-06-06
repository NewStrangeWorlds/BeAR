Setting up a retrieval
======================

To start a retrieval calculation, BeAR requires the following files in the folder the executable is called with:

  - ``retrieval.config`` - the main configuration file for the retrieval
  
  - ``forward_model.config`` - the configuration file of chosen forward forward_model
  
  - ``prior.config`` - the setup list for the prior distributions of the free parameters
  
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
  
  Descriptions of the forward models can be found :ref:`here <sec:opacity_data>` 

| ``Spectral resolution``
|  Sets the constant step in wavenumber space that will be used for the computation of the 
  high-resolution spectrum. Note that this high-resolution grid should generally be finer than that of the
  observational data.

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
  
