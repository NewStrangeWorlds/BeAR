Setting up and starting a retrieval
===================================

A retrieval is started from within its retrieval folder by running the accompanying
Python launcher script:

.. code:: bash

   python run_multinest.py

The launcher script builds the retrieval from the :ref:`configuration files <sec:config_files>` in the
folder and drives the nested-sampling code MultiNest. Internally it uses the ``bear`` pybind module to
construct a ``bear.Config``, ``bear.Retrieval`` and ``bear.PostProcess`` object, converts the MultiNest unit-cube
parameters through the retrieval object, and evaluates the likelihood. After the sampling has finished, the
same script runs the postprocessing step. An example launcher, ``run_multinest.py``, is provided in every
example folder and can be used as a template.

.. code:: python

   #load the retrieval configuration file
   model_config = bear.Config(retrieval_folder)

   #create a BeAR retrieval object
   model = bear.Retrieval(model_config)

   ...

   pymultinest.run(
     loglike,
     priors,
     model.nbParameters(),
     resume = False,
     verbose = True,
     importance_nested_sampling = True,
     sampling_efficiency = 0.8,
     n_live_points = 800,
     max_iter = 0,
     outputfiles_basename = retrieval_folder)

   #create a BeAR retrieval post process object
   post_process = bear.PostProcess(model_config)
   post_process.run()

All sampler settings (importance nested sampling, sampling efficiency, number of live points, maximum number of
iterations, resume, and console feedback) are passed as arguments to ``pymultinest.run`` in the launcher script.
To restart a retrieval that has been interrupted, set ``resume = True`` in the launcher script; all necessary
MultiNest files from the previous run must be present in the retrieval folder. To perform only the postprocessing
step, run a script that creates a ``bear.PostProcess`` object and calls its ``run`` method; the posterior data must
be present in the folder. The GitHub repository of BeAR contains an example for each forward model that can be used
to test the retrieval code or as a template for other retrievals. The example folders contain all necessary files
to run a retrieval calculation.


.. _sec:config_files:

Configuration files
-------------------

BeAR requires the following files in the retrieval folder:

  - ``retrieval.toml`` - the main configuration file for the retrieval

  - ``forward_model.toml`` - the configuration file of the chosen forward model

  - ``priors.config`` - the setup list for the prior distributions of the free parameters

  - ``observations.list`` - the list of observational data files that the retrieval should use

Optionally, the postprocessing step can be configured with the following file:

  - ``post_process.toml`` - the configuration file for the postprocessing step

If this file is not present, BeAR will use default settings for the postprocessing.

The main, forward model, and postprocess configuration files use the `TOML <https://toml.io/>`_ format. The
prior distributions file ``priors.config`` uses the plain-text format described in the
:ref:`section <sec:prior_distributions>` on prior distributions.


Main retrieval file
...................

The file ``retrieval.toml`` contains the basic information for the retrieval setup.

.. include:: ../examples/retrieval.toml_example
   :literal:

The ``[general]`` table contains the following keys:

| ``use_gpu``
|  Boolean that determines the use of the graphics card. Set to ``true`` to run the calculation on the GPU or
  ``false`` to run purely on the CPU.

| ``nb_omp_threads``
|  Integer that sets the number of processor cores used for parallel computing on the CPU. Some parts of BeAR
  will still run on the CPU, even if ``use_gpu`` is enabled. This, for example, is the case for the FastChem
  chemistry code that does not run on GPUs. BeAR will run certain calculations in parallel on the CPU as well,
  using OpenMP. Set this parameter to ``0`` if you want to use all available cores. Note that OpenMP can only
  use a single, multi-core processor.

The ``[retrieval]`` table contains the following keys:

| ``forward_model_type``
|  String that sets the forward model to be used. BeAR currently supports the following models:

     - ``transmission`` - Transmission spectrum

     - ``secondary_eclipse`` - Secondary eclipse / occultation spectrum

     - ``secondary_eclipse_bb`` - Secondary eclipse using a black-body planet

     - ``emission`` - Emission spectrum

     - ``flat_line`` - Fits a flat line to the data

     - ``phase_curve`` - High-resolution phase-curve model

  Short aliases (``trans``, ``se``, ``se_bb``, ``em``, ``fl``, ``pc``) are accepted as well.
  Descriptions of the forward models can be found :ref:`here <sec:forward_models>`.

| ``spectral_discretisation``
|  String that sets the parametrisation of the spectral grid used for the computation of the high-resolution
  spectrum. This high-resolution grid should generally be finer than that of the observational data. The
  following options are available:

     - ``const_wavelength`` - a constant step in wavelength space, with the step size given by
       ``spectral_resolution`` in :math:`\mathrm{\mu m}`

     - ``const_wavenumber`` - a constant step in wavenumber space, with the step size given by
       ``spectral_resolution`` in :math:`\mathrm{cm}^{-1}`

     - ``const_resolution`` - a constant spectral resolution :math:`R = \lambda/\Delta\lambda`, with the value
       given by ``spectral_resolution``

  Note that the key is spelled with the double-s ``discretisation``. In the previous configuration format the
  spectral grid was set with a single line of the form ``keyword value``; this is now split into the two separate
  keys ``spectral_discretisation`` and ``spectral_resolution``.

| ``spectral_resolution``
|  Float that sets the step size or resolution of the spectral grid. Its meaning depends on the chosen
  ``spectral_discretisation`` (see above).

| ``opacity_data_folder``
|  String with the location of the folder that contains the opacities for the gas species. Details on the required
  format of the opacity data can be found in this :ref:`section <sec:opacity_data>`.

| ``use_error_inflation``
|  Boolean that determines the use of the error inflation. This will artificially enlarge the error bars of the
  observational data and, thus, will generally make it easier for the retrieval to find a solution. The use of the
  error inflation acknowledges that certain physical or chemical processes are missing from the simple forward model
  of the retrieval. The form of the employed error inflation is described in
  `Kitzmann et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020ApJ...890..174K/>`_. When this option is set to
  ``true``, BeAR automatically appends an internal error-inflation prior to the parameter list; this prior must **not**
  be added to ``priors.config`` by the user (see the :ref:`section <sec:prior_distributions>` on prior distributions).

The ``[retrieval]`` table also accepts the following optional keys:

| ``spectral_resolution_highres``
|  Float that, when greater than ``0``, enables a separate high-resolution spectral grid for high-resolution
  observations. If it is omitted or set to ``0``, no separate high-resolution grid is created.

| ``forward_model_config`` / ``priors_config`` / ``post_process_config``
|  Strings that override the default file names of the forward model, prior distributions, and postprocess
  configuration files, respectively. If omitted, the defaults ``forward_model.toml``, ``priors.config`` and
  ``post_process.toml`` are used.


Forward model configuration file
................................

The file ``forward_model.toml`` contains the configuration for the forward model.
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

The configuration file for the postprocess step is called ``post_process.toml``. This file is optional and BeAR will use
default settings if it is not present. It contains a set of flat, top-level keys:

.. include:: ../examples/post_process.toml_example
   :literal:

The following keys are available:

| ``delete_sampler_files``
|  Boolean. If ``true``, BeAR deletes the MultiNest files that are not required for the postprocessing step once
  the postprocessing has finished.

| ``save_spectra``
|  Boolean. If ``true``, BeAR computes and writes out the spectra for the posterior sample.

| ``save_temperatures``
|  Boolean. If ``true``, BeAR writes out the temperature structures for the posterior sample.

| ``save_effective_temperatures``
|  Boolean. If ``true``, BeAR computes and writes out the effective temperatures for the posterior sample. This is
  only applicable to emission spectroscopy retrievals.

| ``save_contribution_functions``
|  Boolean. If ``true``, BeAR computes and writes out the contribution functions for the best-fit model.

| ``species_to_save``
|  Array of strings with the chemical species whose mixing ratios should be written out for the posterior sample,
  for example ``["H2O", "CO"]``. An empty array disables this output.

If the ``post_process.toml`` file is absent, BeAR uses default settings for all of these keys.


.. _sec:output_files:

Output files
------------

After a retrieval calculation has been finished, the retrieval folder will contain a set of output files, either directly
from the MultiNest sampler or from the postprocess step of BeAR.

The most important MultiNest files are:

  - ``post_equal_weights.dat`` - the posterior distributions of the model parameters and likelihood values

  - ``summary.dat`` and ``stats.dat`` - a basic summary and some statistics of the nested sampling results, including the
    Bayesian evidence

At the start of a retrieval, BeAR additionally writes the file ``cube_parameters.dat`` to the retrieval folder.
It lists, for each column of the posterior sample, the name of the corresponding parameter and notes which
parameters are fixed (delta priors) or linked to other parameters.

The folder will also contain additional MultiNest files that were used during the nested sampling process. More detailed
descriptions of the files' contents can be found in the MultiNest documentation in its `GitHub repository <https://github.com/farhanferoz/MultiNest/>`_

If the corresponding option in the optional ``post_process.toml`` has been enabled, BeAR will delete MultiNest files that are
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
