Setting up and starting a retrieval
===================================

To start a retrieval calculation, the executable ``bear`` needs to be called from the command line
with a folder that contains the necessary :ref:`configuration files <sec:config_files>` as a command line argument:

.. code:: bash

   ./bear path_to_retrieval_folder/

The folder typically contains all the necessary configuration files for the retrieval and observational data.
:ref:`Output files <sec:output_files>` will also be written to this folder.

Additionally, BeAR can optionally use the following command line arguments:

  - ``-r`` - restart a retrieval that has been interrupted
  
  - ``-p`` - perform only the postprocessing step

The additional arguments are added after the folder path:

.. code:: bash

   ./bear path_to_retrieval_folder/ -p

To restart a retrieval, the folder needs to contain all necessary MultiNest files from the previous run. Likewise, for
the postprocess, the posterior data needs to be present in the folder. The GitHub repository of BeAR contains an example for each 
forward model that can be used to test the retrieval code or as templates for other retrievals. 
The example folders contain all necessary files to run a retrieval calculation.


.. _sec:config_files:

Configuration files
-------------------

BeAR requires the following files in the folder the executable is called with:

  - ``retrieval.config`` - the main configuration file for the retrieval
  
  - ``forward_model.config`` - the configuration file of chosen forward forward_model
  
  - ``priors.config`` - the setup list for the prior distributions of the free parameters
  
  - ``observations.list`` - the list of observational data files that the retrieval should use

Optionally, the postprocessing step can be configured with the following file:

  - ``post_process.config`` - the configuration file for the postprocessing step

If this file is not present, BeAR will use default settings for the postprocessing.


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

The ``priors.config`` file contains the information on the prior distributions of the free parameters. More information
on the format of the prior distributions file can be found in the :ref:`section <sec:prior_distributions>` and in the
description of each forward model.


Observational data file
.......................

The ``observations.list`` file contains a list of data files with the observational data that the retrieval should use.
Its structure depends on the chosen model and is discussed in the :ref:`section <sec:obs_file>` on the observational data.


Postprocess configuration file
..............................

During the postprocess step after a retrieval calculation has been finished, BeAR can perform additional calculations. This
includes the computation of spectra for the posterior sample, writing out all temperature structures, or computing the
effective temperatures for emission spectroscopy retrievals. 

The configuration file for the postprocess step is called
``post_process.config``. This file is optional and BeAR will use default settings if it is not present. The structure of the
file depends on the chosen forward model is discussed in the :ref:`section <sec:forward_models>` on forward models.


.. _sec:output_files:

Output files
------------

After a retrieval calculation has been finished, the retrieval folder will contain a set of output files, either directly
from the MultiNest sampler or from the postprocess step of BeAR.

The most important MultiNest files are:

  - ``post_equal_weights.dat`` - the posterior distributions of the model parameters and likelihood values
  
  - ``summary.dat`` and ``stats.dat`` - a basic summary and some statistics of the nested sampling results, including the
    Bayesian evidence
  
The folder will also contain additional MultiNest files that were used during the nested sampling process. More detailed 
descriptions of the files' contents can be found in the MultiNest documentation in its `GitHub repository <https://github.com/farhanferoz/MultiNest/>`_ 

If the corresponding option in the optional ``post_process.config`` has been enabled, BeAR will delete MultiNest files that are 
not required for the postprocessing step.

The postprocess step will write out additional files, depending on the chosen forward model. This can include:

  - ``spectrum_post_XXXX.dat`` - the spectra for the observation/instrument ``XXXX`` for the posterior sample. Each observational data
    set used in the retrieval will have a separate posterior spectrum file, where  ``XXXX`` is the name stated in the header of the observational 
    data file. The spectra are binned to each observational data. The first column contains the wavelength in :math:`\mathrm{\mu m}`, while all
    other columns are the spectra for each posterior sample. Thus, there are as many spectrum columns as there are posterior samples in the 
    posterior distribution file ``post_equal_weights.dat``. If the original data set was either band-spectroscopy or photometry, the wavelengths
    in the first column refer to the centre of each spectral bin.

  - ``spectrum_best_fit_hr.dat`` - the high-resolution spectrum for the best-fit model, i.e. the model with the highest likelihood. This spectrum
    is saved at the same resulution as the high-resolution grid used in the retrieval. The first column contains the wavelength in :math:`\mathrm{\mu m}`, while
    the second column is the high-resolution spectrum.
  
  - ``temperature_structures.dat`` - the temperature structures for the posterior sample. The first column is the atmospheric pressure in bar. All other columns
    contain the temperatures at these pressures for each posterior sample. The number of temperature columns is equal to the number of posterior samples.
  
  - ``effective_temperatures.dat`` - the effective temperatures for the posterior sample. Each line contains the effective temperature for one posterior sample.

  - ``chem_XXX.dat`` - the mixing ratio of a chemical species ``XXX`` (for example H2O) for the posterior sample. The first column is the atmospheric pressure in bar. 
    All other columns contain the mixing ratios at these pressures for each posterior sample. The number of mixing ratio columns is equal to the number of posterior samples.

  - ``contribution_function_XXXX.dat`` - the contribution functions for the observation/instrument ``XXXX`` for the best-fit model. 
    Each observational data set used in the retrieval will have a separate contribution file, where  ``XXXX`` is the name stated in the header of the observational 
    data file. The first column contains the atmospheric pressure in bar, while all other columns contain the contribution functions for each wavelength/wavelength bin
    of the observational data. The number of contribution function columns is, thus, equal to the number of wavelengths/wavelength bins in the observational data.
    Note that the wavelengths are not saved in this file and have to be taken from either the corresponding observational data file or spectrum posterior file.