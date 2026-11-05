
.. _sec:python_interface_details:

############################################
Detailed description of the Python interface
############################################


Retrieval and PostProcess classes
#################################

These two classes are the main classes of the ``bear`` module. The Retrieval class is used to set up and perform
the retrieval calculations, while the PostProcess class is used to perform post retrieval calculations based
on the posterior samples obtained from a retrieval run.


bear.Retrieval class
********************

This class describes the retrieval object that is used to perform the retrieval calculations. It is mainly responsible for converting
the model parameters into physical values, computing the forward model spectra, and calculating the likelihoods. The class does not
run a sampler itself. Instead, it provides the prior transformation and the likelihood functions of BeAR, which are handed to an
external Python sampler (e.g. MultiNest via PyMultiNest, or dynesty). The class can also return the model spectra for a given set of
parameters that can then be compared to the observational data in Python.

The class has the following structure:

.. code:: python

   class Retrieval
    __init__(
      bear.Config config)
    __init__(
      bear.Config config,
      bear.GenericConfig model_config,
      list[bear.Observation] observations,
      list[bear.Prior] priors)
    bear.SpectralGrid spectral_grid
    int nbParameters()
    pair(list[float], list[float]) convertCubeParameters(
      list[float] cube_parameters)
    list[float] convertToPhysicalParameters(
      list[float] parameters)
    list[float] convertToParameterUnits(
      list[float] physical_parameters)
    float computeLikelihood(
      list[float] physical_parameters)
    bear.ForwardModelOutput computeModel(
      list[float] physical_parameters,
      bool return_high_res_spectrum)
    bear.AtmosphereOutput computeAtmosphereStructure(
      list[float] physical_parameters,
      list[str] species_symbols)
    numpy.ndarray[float32] computeHighResModelFlat(
      list[float] parameters,
      int observation_index = 0,
      bool apply_projection = True)
    bear.HighResLayout highResLayout(
      int observation_index = 0)
    bear.HighResModelExport computeHighResModelMatrices(
      list[float] parameters,
      int observation_index = 0,
      int order_index = 0)

The class contains two different constructors. The first one takes a single argument that is the general configuration
of BeAR:

.. code:: python

   retrieval = bear.Retrieval(config)

In this case, the other configuration files for the chosen forward model, the observational data, and the prior information
are read in during the setup of the object from the folder set in the ``config.retrieval_folder_path`` variable.

The second constructor takes additional arguments for the model configuration, observations, and priors:

.. code:: python

   retrieval = bear.Retrieval(
    config,
    model_config,
    observations,
    priors)

These config objects need to be created beforehand with the corresponding classes of the ``bear`` module as described in the following sections.
The ``model_config`` argument is one of the forward-model configuration classes (e.g. ``bear.TransmissionModelConfig``
or ``bear.EmissionModelConfig``), all of which derive from the common base class ``bear.GenericConfig``.

**Free and full parameter vectors.** BeAR distinguishes between the *free* parameters, i.e. those that are
actually sampled, and the *full* parameter vector that the forward model works with. The full vector contains
one entry for every prior in the order in which the priors are defined, including ``delta`` priors (fixed values)
and ``linked`` priors, which are not sampled. The free vector only contains the sampled parameters. The methods
of the Retrieval class use these two representations as follows:

- ``nbParameters()`` returns the number of *free* parameters.
- ``convertCubeParameters()`` takes a hypercube vector of length ``nbParameters()`` and returns *free* vectors.
- ``convertToParameterUnits()`` takes a *free* vector and returns a *free* vector.
- ``computeLikelihood()``, ``computeHighResModelFlat()``, and ``computeHighResModelMatrices()`` take a *free*
  vector of physical parameters; the delta and linked priors are expanded internally.
- ``convertToPhysicalParameters()`` takes a *free* vector (in prior units) and returns the *full* vector of
  physical parameters.
- ``computeModel()`` and ``computeAtmosphereStructure()`` expect the *full* vector of physical parameters, i.e.
  exactly the output of ``convertToPhysicalParameters()``.

If a retrieval contains neither delta nor linked priors, the free and the full vector are identical, and the
outputs of ``convertCubeParameters()`` can be passed directly to ``computeModel()``. Otherwise, a free vector
of physical parameters can be expanded to the full vector via
``convertToPhysicalParameters(convertToParameterUnits(physical_parameters))``.

The Retrieval class has the following attributes:

``spectral_grid``
  A ``bear.SpectralGrid`` object that contains the spectral grid information for the retrieval.

The Retrieval class has the following methods:

``nbParameters()``
  This method returns an ``int`` that contains the number of *free* parameters of the retrieval, i.e. the number of
  parameters that are actually sampled. Delta and linked priors are not counted. This is the dimension that has to be
  handed to the external sampler.

``convertCubeParameters(cube_parameters)``
  This method takes a ``list[float]`` of hypercube parameters of length ``nbParameters()`` as input and returns a
  ``pair(list[float], list[float])`` that contains the values of the hypercube parameters converted into their prior
  distributions and their physical values. For the latter, the units of the priors are taken into account. Both returned
  lists are *free* vectors of length ``nbParameters()``. This method is typically used as the prior transform in an
  external sampler such as MultiNest.

``convertToPhysicalParameters(parameters)``
  This method takes a ``list[float]`` of *free* parameters in the units of their priors (e.g. one row of a
  ``post_equal_weights.dat`` file without the likelihood column, or the first element returned by
  ``convertCubeParameters()``) and returns the *full* ``list[float]`` of physical parameters, i.e. one value per prior
  including delta and linked priors. The output has the layout expected by ``computeModel()`` and
  ``computeAtmosphereStructure()``.

``convertToParameterUnits(physical_parameters)``
  This method takes a *free* ``list[float]`` of physical parameters (length ``nbParameters()``) and returns a *free*
  ``list[float]`` with the values converted back into the units of the corresponding priors. This is the inverse of the
  physical-value conversion done by ``convertCubeParameters()`` and is used, for example, to write posterior samples in
  the format of the ``post_equal_weights.dat`` file.

``computeLikelihood(physical_parameters)``
  This method takes a *free* ``list[float]`` of physical parameters (the second element returned by
  ``convertCubeParameters()``) as input, expands it internally to the full parameter vector, calculates a forward model
  for these parameters, compares it to the observational data and returns a ``float`` that contains the corresponding
  log-likelihood value for this parameter set.

``computeModel(physical_parameters, return_high_res_spectrum)``
  This method takes the *full* ``list[float]`` of physical parameters (as returned by ``convertToPhysicalParameters()``)
  as input and returns an object of type ``bear.ForwardModelOutput`` that contains the corresponding model spectra.
  If the argument ``return_high_res_spectrum`` is set to ``True``, the high-resolution spectrum will also be returned
  in the output object. Note that, unlike ``computeLikelihood()``, this method does not expand a free parameter vector.

``computeAtmosphereStructure(physical_parameters, species_symbols)``
  This method takes the *full* ``list[float]`` of physical parameters (as returned by ``convertToPhysicalParameters()``)
  as input and returns a ``bear.AtmosphereOutput`` object that contains the corresponding atmosphere structure for these
  parameters. The argument ``species_symbols`` is a ``list[str]`` that contains the chemical symbols of the species that
  should be included in the atmosphere structure output.

``computeHighResModelFlat(parameters, observation_index=0, apply_projection=True)``
  This method is only meaningful for retrievals that include a high-resolution observation (see
  :ref:`here <sec:high_resolution>`) and requires ``config.use_gpu = True``. It takes a *free* ``list[float]`` of
  physical parameters, computes the forward model and returns the high-resolution model for *all* spectral orders of
  the high-resolution observation with index ``observation_index`` as a flat ``numpy.ndarray`` of type ``float32``.
  The model is Doppler-shifted with the velocity parameters :math:`K_\mathrm{p}`, :math:`v_\mathrm{sys}`, and
  :math:`\Delta\phi` contained in the parameter vector. If ``apply_projection`` is ``True``, the same temporal filter
  that is applied to the data in the likelihood is also applied to the model. The layout of the flat array is
  ``[order_offset * nb_exposures + exposure * nb_pixels + pixel]``, where ``order_offset`` is the cumulative number of
  pixels of all preceding orders; the required sizes can be obtained from ``highResLayout()``. An empty array is returned
  if the retrieval has no high-resolution observation, if ``observation_index`` is out of range, or if the GPU is not used.
  This method is intended for building cross-correlation detection maps, where the model is needed for a large number
  of trial velocities.

``highResLayout(observation_index=0)``
  This method returns a ``bear.HighResLayout`` object that describes the high-resolution observation with index
  ``observation_index``: the number of orders and exposures, the number of pixels per order, the orbital phases, the
  wavelengths, and the filtered data in the same flat layout as used by ``computeHighResModelFlat()``. The ``valid``
  attribute of the returned object is ``False`` if there is no such observation.

``computeHighResModelMatrices(parameters, observation_index=0, order_index=0)``
  This method is a diagnostic helper for high-resolution retrievals and requires ``config.use_gpu = True``. It takes a
  *free* ``list[float]`` of physical parameters and returns a ``bear.HighResModelExport`` object that contains, for the
  spectral order ``order_index`` of the high-resolution observation ``observation_index``, the model matrix as used by
  the likelihood both before and after the temporal filter, together with the filtered data it is compared against.
  The ``valid`` attribute of the returned object is ``False`` if there is no high-resolution observation, if one of the
  indices is out of range, or if the GPU is not used.


bear.PostProcess class
**********************

The postprocess class is used to perform post retrieval calculations based on the posterior samples obtained from a retrieval run.
This includes the calculation of model spectra for the posterior samples, effective temperatures, and other derived parameters. 
This class will either read in the posterior sample or take them as an argument.

Although the underlying C++ class is derived from the retrieval class, ``bear.PostProcess`` is exposed in Python as an
*independent* class that is not a subclass of ``bear.Retrieval``. It only re-exposes the subset of the retrieval
methods listed below. In particular, ``nbParameters()``, ``convertCubeParameters()``, ``computeLikelihood()``, and the
high-resolution methods are *not* available on a ``bear.PostProcess`` object.

The class has the following structure:

.. code:: python

    class PostProcess:
      __init__(
        bear.Config config)
      __init__(
        bear.Config config, 
        bear.GenericConfig model_config, 
        bear.GenericConfig post_process_config,
        list[bear.Observation] observations,
        list[bear.Prior] priors)
      bear.SpectralGrid spectral_grid
      list[float] convertToPhysicalParameters(
        list[float] parameters)
      list[float] convertToParameterUnits(
        list[float] physical_parameters)
      bear.ForwardModelOutput computeModel(
        list[float] physical_parameters,
        bool return_high_res_spectrum)
      bear.AtmosphereOutput computeAtmosphereStructure(
        list[float] physical_parameters,
        list[str] species_symbols)
      bool run()
      bool run(
        str posterior_file_path)
      bool run(
        list[list[float]] posteriors,
        list[float] log_likelihoods)

The constructors are similar to those of the ``bear.Retrieval`` class except that the second one takes an
additional argument for the postprocess configuration. Both ``model_config`` and ``post_process_config`` are
objects of the forward-model specific configuration classes described below (e.g. ``bear.EmissionModelConfig`` and
``bear.EmissionPostProcessConfig``), which all derive from ``bear.GenericConfig``.

The attribute ``spectral_grid`` and the methods ``convertToPhysicalParameters()``, ``convertToParameterUnits()``,
``computeModel()``, and ``computeAtmosphereStructure()`` behave exactly as described for the ``bear.Retrieval`` class,
including the distinction between free and full parameter vectors: ``convertToPhysicalParameters()`` converts a free
vector in prior units (e.g. one row of a ``post_equal_weights.dat`` file) into the full vector of physical parameters,
which is the input expected by ``computeModel()`` and ``computeAtmosphereStructure()``.

The PostProcess class has the following additional methods, all of which return ``True`` on success and ``False``
if an error occurred:

``run()``
  This method starts the post process calculations using the standard MultiNest output file ``post_equal_weights.dat``
  located in the folder set in the ``config.output_path`` variable. 
  The output will be saved in the folder set by the ``config.post_output_path`` variable.

``run(posterior_file_path)``
  This method takes a ``str`` argument that is the path to a file with the posterior samples that should be
  post-processed. The file has to have the format of the MultiNest ``post_equal_weights.dat`` file, i.e. one row per
  sample containing the free parameters in the units of their priors followed by the log-likelihood as the last column.
  The output will be saved in the folder set by the ``config.post_output_path`` variable.

``run(posteriors, log_likelihoods)``
  This method takes a ``list[list[float]]`` of posterior samples and a ``list[float]`` of corresponding log-likelihoods
  as arguments and starts the post process calculations based on these samples. Each sample is a *free* parameter
  vector of length equal to the number of free parameters in the units of the priors, i.e. the same format as the
  rows of the ``post_equal_weights.dat`` file; the conversion to physical parameters is done internally. This overload
  is used when the posterior is produced by a Python sampler (e.g. dynesty or nautilus) without a MultiNest output file.
  The output will be saved in the folder set by the ``config.post_output_path`` variable.


bear.Config
###########

This class describes the general model configuration of BeAR and is equivalent to the information located in a
``retrieval.toml`` file. The class is used to load the configuration files of BeAR and to set up the retrieval object.

It has the following attributes and methods:

.. code:: python

   class Config:
       __init__()
       __init__(str file_path)
       __init__(
          bool use_gpu,
          str forward_model_type,
          str cross_section_file_path,
          str spectral_discretisation,
          float resolution,
          str output_folder,
          str post_process_output_folder)
      bool loadConfigFile(str retrieval_folder)
      forward_model_type : str = ""
      retrieval_folder_path : str = ""
      output_path : str = ""
      post_output_path : str = ""
      wavenumber_file_path : str = ""
      cross_section_file_path : str = ""
      spectral_disecretisation : int = 0
      spectral_resolution : float = 0
      use_error_inflation : bool = False
      use_gpu : bool = False
      nb_omp_processes : int = 0

The class contains three different constructors. The first one is the default constructor that does not take any arguments.

.. code:: python

   model_config = bear.Config()

In this case the object will be created will all values set to their default values. The second constructor takes a string
argument that is the path to the retrieval folder. This is the folder where the configuration file ``retrieval.toml`` is
located.

.. code:: python

   model_config = bear.Config(folder)

The function will load the configuration file and set the attributes of the object accordingly. The third constructor is the
most flexible one and allows the user to set the most important attributes of the object directly. The arguments are:

.. code:: python

   model_config = bear.Config(
    bool use_gpu,
    str forward_model_type,
    str cross_section_file_path,
    str spectral_discretisation,
    float spectral_resolution,
    str output_folder,
    str post_process_output_folder)

All other attributes are set to their default values. The attributes of the object can be accessed and changed directly after
the object has been created. For example, the forward model type can be changed via:

.. code:: python

   model_config.forward_model_type = "emission"

or the error inflation can be enabled via:

.. code:: python

   model_config.use_error_inflation = True

The Config class has the following attributes:

``forward_model_type``
  A ``str`` that contains the choice of the forward model. The available model types can be found :ref:`here <sec:forward_models>`

``retrieval_folder_path``
  A ``str`` that contains the path to the folder where the configuration files are located.

``output_path``
  A ``str`` that contains the path to the folder where the sampler (e.g. MultiNest) output files will be stored.

``post_output_path``
  A ``str`` that contains the path to the folder where the post process output files will be stored.

``wavenumber_file_path``
  A ``str`` that contains the path to the file that contains the wavenumber grid. If the standard HELIOS-k opacities 
  with a constant wavenumber step of 0.01 cm^-1 are used, this file is not needed.

``cross_section_file_path``
  A ``str`` that contains the path to the file that contains the opacity data.

``spectral_disecretisation``
  An ``int`` that contains the choice of the spectral discretisation: ``0`` for ``const_wavenumber``, ``1`` for
  ``const_wavelength``, and ``2`` for ``const_resolution``. The meaning of these options is described
  :ref:`here <sec:config_files>`. Only the seven-argument constructor accepts the string form; when the attribute is set
  directly, the integer code has to be used. Note that the attribute name is misspelled in the wrapper
  (``disecretisation``) and has to be used with exactly this spelling.

``spectral_resolution``
  A ``float`` that contains the spectral resolution or constant wavelength/wavenumber step of the data.

``use_error_inflation``
  A ``bool`` that determines if the error inflation is enabled.

``use_gpu``
  A ``bool`` that determines if the GPU is used for the calculations.

``nb_omp_processes``
  An ``int`` that contains the number of OpenMP processes that are used for the calculations.
  A number of 0 implies that the maximum available number of processes will be used.

The Config class has the following method:

``loadConfigFile(retrieval_folder)``
  This method takes a ``str`` argument that is the path to the retrieval folder containing the ``retrieval.toml`` file,
  reads the file and sets the attributes of the object accordingly. It returns a ``bool`` that is ``True`` on success.
  The second constructor calls this method internally.

.. note::

   The underlying configuration also contains further options that are read from ``retrieval.toml`` but are
   not exposed as writable attributes of the Python ``Config`` class. These include ``spectral_resolution_highres``
   (the spectral resolution of the high-resolution grid, used by the phase-curve forward model) and
   ``nb_disort_streams`` (the number of streams used by the DisORT++ radiative transfer solver). Set these
   through the ``retrieval.toml`` file.


bear.SpectralGrid
#################

This class describes the spectral grid that is used for the retrieval calculations. The class has the following structure:

.. code:: python
  
   class SpectralGrid
    __init__(
      bear.Config config)
    __init__(
      bear.Config config,
      float wavelength_min,
      float wavelength_max)
    wavenumber_list : list[float]
    wavelength_list : list[float]
    float wavelengthToWavenumber(
      float wavelength)

The class contains two constructors. The first one only takes the general config as an argument:

.. code:: python

   spectral_grid = bear.SpectralGrid(config)

In this case, only the full wavenumber grid of the opacity data (from ``config.wavenumber_file_path`` or the standard
HELIOS-K grid) is loaded; the actual spectral grid used for the model calculations is left empty and is set up later
from the wavelength ranges of the observations. The second constructor additionally takes the minimum and maximum
wavelength of the spectral grid in units of micron:

.. code:: python

   spectral_grid = bear.SpectralGrid(
    config,
    wavelength_min,
    wavelength_max)

The spectral grid will be created based on the provided wavelength range and the spectral discretisation and resolution
set in the configuration. This is the constructor to use when calling a forward model directly (see the section on
the direct forward-model bindings at the end of this page). Note that the retrieval and post process classes will
automatically create the spectral grid during their setup, so the user usually does not need to create this object directly.

The SpectralGrid class has the following attributes and methods:

``wavenumber_list``
  A ``list[float]`` that contains the wavenumber grid in units of cm^-1.

``wavelength_list``
  A ``list[float]`` that contains the wavelength grid in units of micron.

``float wavelengthToWavenumber(float wavelength)``
  This method takes a ``float`` argument that is a wavelength in units of micron and returns a ``float``
  that contains the corresponding wavenumber in units of cm^-1.


bear.Observation
################

This class describes the a single observation that are usually read in from a file with a structure described :ref:`here <sec:observational_data>`.
The class has the following structure:

.. code:: python

   class Observation:
    __init__()
    __init__(
      str name,
      str type)
    __init__(
      str name, 
      str type, 
      list[float] wavelengths,
      list[float] observation_data,
      list[float] observation_error)
    __init__(
      str name, 
      str type, 
      list[list[float]] wavelength_bins,
      list[float] observation_data,
      list[float] observation_error)
    name : str = ""
    type : str = ""
    wavelengths : list[float]
    bin_wavelength_edges : list[list[float]]
    data : list[float]
    data_error : list[float]
    likelihood_weight : list[float]
    instrument_profile_fwhm : list[float]
    filter_response : list[list[float]]
    filter_detector_type : str = ""
    spectrum_modifier_id : str = ""

The class contains four different constructors. The first one is the default constructor that does not take any arguments.

.. code:: python

   observation = bear.Observation()

In this case the object will be created will all values set to their default values. The second constructor takes two string
argumentse that is the name of the observation and its type:

.. code:: python

   observation = bear.Observation(name, type)

The third constructor is used for spectroscopy observations that are given as a list of wavelengths, observation data, 
and observation errors:

.. code:: python

   observation = bear.Observation(
    name, 
    type, 
    wavelengths, observation_data, 
    observation_error)

Finally, the fourth constructor is used for photometry or band-spectroscopy observations that are given as a list of wavelength bins, observation data, and observation errors:

.. code:: python

   observation = bear.Observation(
    name, 
    type, 
    wavelength_bins, 
    observation_data, 
    observation_error)

The variable ``wavelength_bins`` is a two-dimensional array that contains the wavelengths of the  bin edges.

All attributes can also be set directly after the object has been created. For example, the observation type can be changed via:

.. code:: python

   observation.type = "photometry"

The Observation class has the following attributes:

``name``
  A ``str`` that contains the name of the observation.

``type``
  A ``str`` that contains the type of the observation. The available options can be found :ref:`here <sec:observational_data>`

``wavelengths``
  A ``list[float]`` that contains the wavelengths of the observation. This value needs to be set if the observation is of type ``spectroscopy``. 
  For all other types, this value is ignored.

``bin_wavelength_edges``
  A ``list[list[float]]`` that contains the wavelength edges of the observation bins. This value needs to be set if the observation is of 
  type ``photometry`` or ``band-spectroscopy``. It is ignored for all other types.

``data``
  A ``list[float]`` that contains the observation data.

``data_error``
  A ``list[float]`` that contains the observation errors.

``likelihood_weight``
  A ``list[float]`` that contains the likelihood weights of the observation. If no weights are provided, they are set to unity.

``instrument_profile_fwhm``
  A ``list[float]`` that contains the optional full width at half maximum of the instrument profile.

``filter_response``
  A ``list[list[float]]`` that contains the optional filter response of the observation. This is a two-dimensional array where the first
  column contains the wavelengths and the second column the response function. The response function will internally be interpolated to the
  wavenumber grid of the retrieval calculation.

``filter_detector_type``
  A ``str`` that contains the type of the detector if a filter function is used.

``spectrum_modifier_id``
  A ``str`` that contains the optional spectrum modifier.


bear.Prior
##########

This class describes a specific model prior distribution used in the retrieval as described :ref:`here <sec:prior_distributions>`.
The class has the following structure:

.. code:: python

   class Prior
    __init__(
      str type, 
      str description, 
      list[float] parameter)
    __init__(
      str type, 
      str description, 
      list[float] parameter,
      str unit)
    type : str = ""
    description : str = ""
    parameter : list[float]
    unit : str = ""
    link_target : str = ""

Note that the *type* of the prior is the first constructor argument and its *description* (the parameter name) the
second one, i.e. the same order as the columns in the ``priors.config`` file. The class contains two different
constructors. The first one is used for priors that do not have a unit:

.. code:: python

   prior = bear.Prior(
    type, 
    description, 
    parameter)

For example, a uniform prior for the surface gravity can be created via:

.. code:: python

   prior = bear.Prior(
    "uniform", 
    "log_g", 
    [2.0, 5.0]) 

and a Gaussian prior via:

.. code:: python

   prior = bear.Prior(
    "gaussian", 
    "log_g", 
    [2.99, 0.03]) 

The second constructor is used for priors that have a unit:

.. code:: python

   prior = bear.Prior(
    type, 
    description, 
    parameter,
    unit)

For example, a Gaussian prior in solar radii can be created via:

.. code:: python

   prior = bear.Prior(
    "gaussian", 
    "r_star", 
    [1.57, 0.07],
    "Rs")

All attributes can also be set directly after the object has been created. For example, the unit can be changed via:

.. code:: python

   prior.unit = "Rj"

The Prior class has the following attributes:

``type``
  A ``str`` that contains the type of the prior. The available options can be found :ref:`here <sec:prior_distributions>`

``description``
  A ``str`` that contains the description/name of the prior. This is the parameter name that is matched against the
  parameter names of the forward model.

``parameter``
  A ``list[float]`` that contains the parameters of the prior. For a prior that only requires one parameter, 
  this variable still needs to be a list.

``unit``
  A ``str`` that contains the optional unit of the prior. A list of possible units can be found :ref:`here <sec:prior_distributions>`

``link_target``
  A ``str`` that contains the name (``description``) of another prior to which this prior is linked. It is only used
  for priors of type ``linked`` and ties the value of this parameter to that of the target parameter instead of
  sampling it independently. There is no constructor argument for this attribute; it has to be set directly after
  the object has been created, e.g. ``prior.link_target = "temperature_1"``.



bear.ForwardModelOutput
#######################

  .. code:: python

   class ForwardModelOutput:
      neglect_model : bool
      spectrum_obs : list[list[float]]
      spectrum : list[float]

This structure is used to store the output of a forward model calculation. It has the following attributes:

``neglect_model``
  A ``bool`` that determines if the model should be neglected in the likelihood calculations. This is set to ``True``
  if the forward model calculation failed. This can occur, for example, if unphysical parameter values were used, such as when the sum of all
  volume mixing ratios exceeds unity.

``spectrum_obs``
  A ``list[list[float]]`` that contains the binned spectrum for the observational data. 
  The first index refers to the observation and the second one to the wavelengths of each observation. 
  The size of the first dimension corresponds to the number of observations used in the retrieval. 
  The wavelength points for the second dimension are identical to those read in from the observation files and are not stored here but can be accessed
  via the ``bear.Observation`` objects.

``spectrum``
  A ``list[float]`` that contains the high-resolution spectrum calculated on the spectral grid of the retrieval.


bear.AtmosphereOutput
#####################

  .. code:: python

   class AtmosphereOutput:
      neglect_model : bool
      pressure : list[float]
      altitude : list[float]
      temperature : list[float]
      species_symbols : list[str]
      mixing_ratios : list[list[float]]

This structure is used to store the output of an atmosphere structure calculation. It has the following attributes:

``neglect_model``
  A ``bool`` that determines if the model should be neglected in the likelihood calculations. This is set to ``True``
  if the atmosphere structure calculation failed. This can occur, for example, if unphysical parameter values were used, such as when the sum of all
  volume mixing ratios exceeds unity.

``pressure``
  A ``list[float]`` that contains the pressure values in bar for each level of the atmosphere.

``altitude``
  A ``list[float]`` that contains the altitude values in units of cm for each level of the atmosphere.

``temperature``
  A ``list[float]`` that contains the temperature values for each level of the atmosphere.

``species_symbols``
  A ``list[str]`` that contains the chemical symbols of the species present in the output.

``mixing_ratios``
  A ``list[list[float]]`` that contains the mixing ratios of chemical species for each level of the atmosphere. The first dimension refers to the species
  and the second dimension to the atmospheric levels. The order of the species corresponds to that given in the ``species_symbols`` list.


bear.HighResLayout
##################

  .. code:: python

   class HighResLayout:
      valid : bool
      nb_orders : int
      nb_exposures : int
      kp_ref : float
      vsys_ref : float
      order_nb_pixels : list[int]
      orbital_phases : list[float]
      wavelengths : list[float]
      data_filtered : numpy.ndarray[float32]

This structure describes the layout of a :ref:`high-resolution observation <sec:highres_observations>` and is returned
by ``Retrieval.highResLayout()``. It contains everything that does not change with the trial velocities when building,
for example, a cross-correlation detection map, and is meant to be fetched once and then paired with the output of
``Retrieval.computeHighResModelFlat()``. All attributes are read-only. The structure has the following attributes:

``valid``
  A ``bool`` that is ``False`` if the retrieval has no high-resolution observation with the requested index. In this
  case, all other attributes are empty or zero.

``nb_orders``
  An ``int`` with the number of spectral orders of the observation.

``nb_exposures``
  An ``int`` with the number of exposures (time steps) of the observation.

``kp_ref``
  A ``float`` with the reference planetary orbital velocity :math:`K_\mathrm{p}` of the observation in km/s.
  The retrieved :math:`K_\mathrm{p}` parameter is an offset with respect to this value.

``vsys_ref``
  A ``float`` with the reference systemic velocity :math:`v_\mathrm{sys}` of the observation in km/s.
  The retrieved :math:`v_\mathrm{sys}` parameter is an offset with respect to this value.

``order_nb_pixels``
  A ``list[int]`` with the number of wavelength pixels of each spectral order.

``orbital_phases``
  A ``list[float]`` with the orbital phase of each exposure.

``wavelengths``
  A ``list[float]`` with the wavelengths of all orders in nm, concatenated order by order, i.e. with the layout
  ``[order_offset + pixel]`` where ``order_offset`` is the cumulative number of pixels of all preceding orders.

``data_filtered``
  A one-dimensional ``numpy.ndarray`` of type ``float32`` with the filtered observational data of all orders and
  exposures, in the flat layout ``[order_offset * nb_exposures + exposure * nb_pixels + pixel]``. This is the same
  layout as the array returned by ``Retrieval.computeHighResModelFlat()``, so that the two can be compared element by
  element.


bear.HighResModelExport
#######################

  .. code:: python

   class HighResModelExport:
      valid : bool
      nb_exposures : int
      nb_pixels : int
      wavelengths : list[float]
      model_unfiltered : list[list[float]]
      model_filtered : list[list[float]]
      data_filtered : list[list[float]]

This structure is returned by ``Retrieval.computeHighResModelMatrices()`` and contains, for a single spectral order
of a high-resolution observation, the model as it is built by the likelihood calculation before and after the
temporal filter, alongside the filtered data it is compared against. It is intended as a diagnostic tool; nothing in
it feeds back into the sampling. All attributes are read-only. The structure has the following attributes:

``valid``
  A ``bool`` that is ``False`` if the retrieval has no high-resolution observation, if the observation or order index
  is out of range, or if the GPU is not used. In this case, all other attributes are empty or zero.

``nb_exposures``
  An ``int`` with the number of exposures of the observation.

``nb_pixels``
  An ``int`` with the number of wavelength pixels of the requested order.

``wavelengths``
  A ``list[float]`` with the wavelengths of the requested order in nm.

``model_unfiltered``
  A ``list[list[float]]`` with the Doppler-shifted model before the temporal filter is applied. The first index refers
  to the exposure and the second one to the pixel.

``model_filtered``
  A ``list[list[float]]`` with the model after the temporal filter has been applied, i.e. the quantity that enters the
  likelihood. The indices are the same as for ``model_unfiltered``.

``data_filtered``
  A ``list[list[float]]`` with the filtered observational data of the requested order, with the same indices as
  ``model_filtered``.


Forward Model Config and PostProcess Config
###########################################

Every forward model in BeAR has its own configuration that is used to set up the model and an additional one
for the post process calculations. The configuration parameters are different for every model type and are
described in the :ref:`section <sec:forward_models>` on forward models. Below, the Python classes for the
transmission, secondary-eclipse, secondary-eclipse-blackbody, emission, and flat-line forward models are
described. The high-resolution phase-curve forward model is also available (via the ``phase_curve``
forward-model type), but does not currently expose a dedicated Python configuration class; it is set up
through its ``forward_model.toml`` and ``post_process.toml`` files.

Transmission Spectrum
*********************


Model Config
============

The configuration class for the transmission model has the following structure:

.. code:: python

   class TransmissionModelConfig:
    __init__(str file_path)
    __init__(
      int nb_grid_points,
      float atmos_bottom_pressure,
      float atmos_top_pressure,
      str temperature_profile_model,
      list[str] temperature_profile_parameters,
      list[str] chemistry_model,
      list[list[str]] chemistry_parameters,
      list[str] opacity_species_symbol,
      list[str] opacity_species_folder)
    __init__(
      bool fit_mean_molecular_weight,
      bool fit_scale_height,
      bool use_variable_gravity,
      int nb_grid_points,
      float atmos_bottom_pressure,
      float atmos_top_pressure,
      str temperature_profile_model,
      list[str] temperature_profile_parameters,
      list[str] chemistry_model,
      list[list[str]] chemistry_parameters,
      list[str] opacity_species_symbol,
      list[str] opacity_species_folder,
      list[str] cloud_model,
      list[list[str]] cloud_model_parameters,
      list[str] modules,
      list[list[str]] modules_parameters)
    nb_grid_points : int = 0
    atmos_boundaries : list[float]
    fit_mean_molecular_weight : bool = False
    fit_scale_height : bool = False
    temperature_profile_model : str = ""
    temperature_profile_parameters : list[str]
    chemistry_model : list[str]
    chemistry_parameters : list[list[str]]
    cloud_model : list[str]
    cloud_model_parameters : list[list[str]]
    modules : list[str]
    modules_parameters : list[list[str]]
    opacity_species_symbol : list[str]
    opacity_species_folder : list[str]

The class contains three different constructors. The first one takes a string argument that is the path to the
transmission model configuration file: `file_path`. 

.. code:: python

   model_config = bear.TransmissionModelConfig(file_path)

The second constructor takes the most important parameters of the model as arguments:

.. code:: python

   model_config = bear.TransmissionModelConfig(
      int nb_grid_points,
      float atmos_bottom_pressure,
      float atmos_top_pressure,
      str temperature_profile_model,
      list[str] temperature_profile_parameters,
      list[str] chemistry_model,
      list[list[str]] chemistry_parameters,
      list[str] opacity_species_symbol,
      list[str] opacity_species_folder)

The constructor has the following arguments:

``nb_grid_points``
  An ``int`` that contains the number of grid points of the atmosphere.

``atmos_bottom_pressure``
  A ``float`` that contains the bottom pressure of the atmosphere in bar.

``atmos_top_pressure``
  A ``float`` that contains the top pressure of the atmosphere in bar.

``temperature_profile_model``
  A ``str`` that contains the choice of the temperature profile model. The available options can be found :ref:`here <sec:temperature_profiles>`.

``temperature_profile_parameters``
  A ``list[str]`` that contains the parameters of the temperature profile model.

``chemistry_model``
  A ``list[str]`` that contains the choice of the chemistry model(s). The available options can be found :ref:`here <sec:chemistry_models>`.

``chemistry_parameters``
  A ``list[list[str]]`` that contains the parameters of the chemistry model(s).

``opacity_species_symbol``
  A ``list[str]`` that contains the symbols/formulas of the opacity species used in the radiative transfer calculations.

  .. warning::

     In the current Python wrapper, the *attribute* ``opacity_species_symbol`` of ``bear.TransmissionModelConfig`` is
     erroneously bound to the ``modules`` member of the underlying C++ class. Reading or writing
     ``model_config.opacity_species_symbol`` therefore accesses the list of modules instead of the opacity species.
     The constructor argument is not affected; set the opacity species through the constructor and do not modify this
     attribute afterwards.

``opacity_species_folder``
  A ``list[str]`` that contains the folder where opacity for each species are stored.

The third constructor is the most flexible one and takes additional arguments for advanced model options:

.. code:: python

   model_config = bear.TransmissionModelConfig(
      bool fit_mean_molecular_weight,
      bool fit_scale_height,
      bool use_variable_gravity,
      int nb_grid_points,
      float atmos_bottom_pressure,
      float atmos_top_pressure,
      str temperature_profile_model,
      list[str] temperature_profile_parameters,
      list[str] chemistry_model,
      list[list[str]] chemistry_parameters,
      list[str] opacity_species_symbol,
      list[str] opacity_species_folder,
      list[str] cloud_model,
      list[list[str]] cloud_model_parameters,
      list[str] modules,
      list[list[str]] modules_parameters)

The additional arguments are:

``fit_mean_molecular_weight``
  A ``bool`` that determines if the mean molecular weight is fitted as an additional free parameter.

``fit_scale_height``
  A ``bool`` that determines if the scale height is fitted as an additional free parameter.

``use_variable_gravity``
  A ``bool`` that determines if the variable gravity with altitude is used in the calculations.

``cloud_model``
  A ``list[str]`` that contains the choice of the cloud model(s). The available options can be found :ref:`here <sec:cloud_models>`.

``cloud_model_parameters``
  A ``list[list[str]]`` that contains the parameters of the cloud model(s).

``modules``
  A ``list[str]`` that contains the choice of additional model modules. The available options can be found :ref:`here <sec:forward_models>`.

``modules_parameters``
  A ``list[list[str]]`` that contains the parameters of the additional model module(s).

All attributes can also be set directly after the object has been created. For example, the number of grid points can be changed via:

.. code:: python

   model_config.nb_grid_points = 128


PostProcess Config
==================

The configuration class for the transmission model post process has the following structure:

.. code:: python

   class TransmissionPostProcessConfig:
    __init__(str file_path)
    __init__(
      bool save_temperatures, 
      bool save_spectra, 
      list[str] chemical_species_symbols)
    save_temperatures : bool = False
    save_spectra : bool = True
    delete_sampler_files : bool = False

The class contains two different constructors. The first one takes a string argument that is the path to the
transmission model post process configuration file: `file_path`.

.. code:: python

   post_process_config = bear.TransmissionPostProcessConfig(file_path)

The second constructor takes the most important parameters of the post process as arguments:

.. code:: python

   post_process_config = bear.TransmissionPostProcessConfig(
      bool save_temperatures,
      bool save_spectra,
      list[str] chemical_species_symbols)

The different options are discussed :ref:`here <sec:forward_model_transmission>`.


Secondary Eclipse Spectrum
**************************

Model Config
============

The configuration class for the secondary eclipse model has the following structure:

.. code:: python

   class OccultationConfig:
    __init__(str file_path)
    __init__(
      int nb_grid_points,
      float atmos_bottom_pressure,
      float atmos_top_pressure,
      str temperature_profile_model,
      list[str] temperature_profile_parameters,
      str radiative_transfer_model,
      list[str] radiative_transfer_parameters,
      list[str] chemistry_model,
      list[list[str]] chemistry_parameters,
      list[str] opacity_species_symbol,
      list[str] opacity_species_folder,
      str stellar_spectrum_model,
      list[str] stellar_model_parameters,
      list[str] cloud_model,
      list[list[str]] cloud_model_parameters)
    nb_grid_points : int = 0
    atmos_boundaries : list[float]
    temperature_profile_model : str = ""
    temperature_profile_parameters : list[str]
    radiative_transfer_model : str = ""
    radiative_transfer_parameters : list[str]
    stellar_spectrum_model : str = ""
    stellar_model_parameters : list[str]
    chemistry_model : list[str] = []
    chemistry_parameters : list[list[str]]
    cloud_model : list[str] = []
    cloud_model_parameters : list[list[str]]
    opacity_species_symbol : list[str]
    opacity_species_folder : list[str]

The class contains two different constructors. The first one takes a string argument that is the path to the
secondary eclipse model configuration file: `file_path`.

.. code:: python

   model_config = bear.OccultationConfig(file_path)

The second constructor takes the most important parameters of the model as arguments:

.. code:: python

   model_config = bear.OccultationConfig(
      int nb_grid_points,
      float atmos_bottom_pressure,
      float atmos_top_pressure,
      str temperature_profile_model,
      list[str] temperature_profile_parameters,
      str radiative_transfer_model,
      list[str] radiative_transfer_parameters,
      list[str] chemistry_model,
      list[list[str]] chemistry_parameters,
      list[str] opacity_species_symbol,
      list[str] opacity_species_folder,
      str stellar_spectrum_model,
      list[str] stellar_model_parameters,
      list[str] cloud_model,
      list[list[str]] cloud_model_parameters)

The constructor has the following arguments:

``nb_grid_points``
  An ``int`` that contains the number of grid points of the atmosphere.

``atmos_bottom_pressure``
  A ``float`` that contains the bottom pressure of the atmosphere in bar.

``atmos_top_pressure``
  A ``float`` that contains the top pressure of the atmosphere in bar.

``temperature_profile_model``
  A ``str`` that contains the choice of the temperature profile model. The available options can be found :ref:`here <sec:temperature_profiles>`.

``temperature_profile_parameters``
  A ``list[str]`` that contains the parameters of the temperature profile model.

``radiative_transfer_model``
  A ``str`` that contains the choice of the radiative transfer model. The available options can be found :ref:`here <sec:forward_model_se>`.

``radiative_transfer_parameters``
  A ``list[str]`` that contains the parameters of the radiative transfer model.

``chemistry_model``
  A ``list[str]`` that contains the choice of the chemistry model(s). The available options can be found :ref:`here <sec:chemistry_models>`.

``chemistry_parameters``
  A ``list[list[str]]`` that contains the parameters of the chemistry model(s).

``opacity_species_symbol``
  A ``list[str]`` that contains the symbols/formulas of the opacity species used in the radiative transfer calculations.

``opacity_species_folder``
  A ``list[str]`` that contains the folder where opacity for each species are stored.

``stellar_spectrum_model``
  A ``str`` that contains the choice of the stellar spectrum model. The available options can be found :ref:`here <sec:stellar_spectra>`.

``stellar_model_parameters``
  A ``list[str]`` that contains the parameters of the stellar spectrum model.

``cloud_model``
  A ``list[str]`` that contains the choice of the cloud model(s). The available options can be found :ref:`here <sec:cloud_models>`.

``cloud_model_parameters``
  A ``list[list[str]]`` that contains the parameters of the cloud model(s).


PostProcess Config
==================

The configuration class for the secondary eclipse model post process has the following structure:

.. code:: python

   class OccultationPostProcessConfig:
    __init__(str file_path)
    __init__(
      bool save_temperatures, 
      bool save_spectra, 
      bool save_contribution_functions,
      list[str] chemical_species_symbols)
    save_temperatures : bool = True
    save_spectra : bool = True
    save_contribution_functions : bool = False,
    delete_sampler_files : bool = False

The class contains two different constructors. The first one takes a string argument that is the path to the
secondary eclipse model post process configuration file: `file_path`.

.. code:: python

   post_process_config = bear.OccultationPostProcessConfig(file_path)

The second constructor takes the most important parameters of the post process as arguments:

.. code:: python

   post_process_config = bear.OccultationPostProcessConfig(
      bool save_temperatures,
      bool save_spectra,
      bool save_contribution_functions,
      list[str] chemical_species_symbols)

The different options are discussed :ref:`here <sec:forward_model_se>`.


Secondary Eclipse Spectrum with Planetary Blackbody
***************************************************

Model Config
============

The configuration class for this special secondary eclipse model has the following structure:

.. code:: python

   class OccultationBlackBodyConfig:
    __init__(str file_path)
    __init__(
      str stellar_spectrum_model,
      list[str] stellar_model_parameters)
    stellar_spectrum_model : str = ""
    stellar_model_parameters : list[str]

The class contains two different constructors. The first one takes a string argument that is the path to the
secondary eclipse blackbody model configuration file: `file_path`. 

.. code:: python

   model_config = bear.OccultationBlackBodyConfig(file_path)

The second constructor takes the most important parameters of the model as arguments:

.. code:: python

   model_config = bear.OccultationBlackBodyConfig(
      str stellar_spectrum_model,
      list[str] stellar_model_parameters)

The constructor has the following arguments:

``stellar_spectrum_model``
  A ``str`` that contains the choice of the stellar spectrum model. The available options can be found :ref:`here <sec:stellar_spectra>`.

``stellar_model_parameters``
  A ``list[str]`` that contains the parameters of the stellar spectrum model.

PostProcess Config
==================

The configuration class for the secondary eclipse model post process has the following structure:

.. code:: python

   class OccultationBlackBodyPostConfig:
    __init__(str file_path)
    __init__(bool save_spectra)
    save_spectra : bool = True
    delete_sampler_files : bool = False

The class contains two different constructors. The first one takes a string argument that is the path to the
secondary eclipse model post process configuration file: `file_path`.

.. code:: python

   post_process_config = bear.OccultationBlackBodyPostConfig(file_path)

The second constructor takes the most important parameter of the post process as an argument:

.. code:: python

   post_process_config = bear.OccultationBlackBodyPostConfig(
      bool save_spectra)

Note that, unlike the other post process configuration classes, the name of this class does not contain
``Process``; it is ``bear.OccultationBlackBodyPostConfig``.

The different options are discussed :ref:`here <sec:forward_model_se_bb>`.


Emission Spectrum
*****************

Model Config
============

The configuration class for the emission model has the following structure:

.. code:: python

   class EmissionModelConfig:
    __init__(str file_path)
    __init__(
      int nb_grid_points,
      float atmos_bottom_pressure,
      float atmos_top_pressure,
      str temperature_profile_model,
      list[str] temperature_profile_parameters,
      str radiative_transfer_model,
      list[str] radiative_transfer_parameters,
      list[str] chemistry_model,
      list[list[str]] chemistry_parameters,
      list[str] opacity_species_symbol,
      list[str] opacity_species_folder)
    __init__(
      int nb_grid_points,
      float atmos_bottom_pressure,
      float atmos_top_pressure,
      str temperature_profile_model,
      list[str] temperature_profile_parameters,
      str radiative_transfer_model,
      list[str] radiative_transfer_parameters,
      list[str] chemistry_model,
      list[list[str]] chemistry_parameters,
      list[str] opacity_species_symbol,
      list[str] opacity_species_folder,
      list[str] cloud_model,
      list[list[str]] cloud_model_parameters)
    nb_grid_points : int = 0
    atmos_boundaries : list[float]
    temperature_profile_model : str = ""
    temperature_profile_parameters : list[str]
    radiative_transfer_model : str = ""
    radiative_transfer_parameters : list[str]
    chemistry_model : list[str] = []
    chemistry_parameters : list[list[str]]
    cloud_model : list[str] = []
    cloud_model_parameters : list[list[str]]
    opacity_species_symbol : list[str]
    opacity_species_folder : list[str]

The class contains three different constructors. The first one takes a string argument that is the path to the
emission model configuration file: `file_path`.

.. code:: python

   model_config = bear.EmissionModelConfig(file_path)

The second constructor takes the most important parameters of the model as arguments:

.. code:: python

   model_config = bear.EmissionModelConfig(
      int nb_grid_points,
      float atmos_bottom_pressure,
      float atmos_top_pressure,
      str temperature_profile_model,
      list[str] temperature_profile_parameters,
      str radiative_transfer_model,
      list[str] radiative_transfer_parameters,
      list[str] chemistry_model,
      list[list[str]] chemistry_parameters,
      list[str] opacity_species_symbol,
      list[str] opacity_species_folder,
      list[str] cloud_model,
      list[list[str]] cloud_model_parameters)

The constructor has the following arguments:

``nb_grid_points``
  An ``int`` that contains the number of grid points of the atmosphere.

``atmos_bottom_pressure``
  A ``float`` that contains the bottom pressure of the atmosphere in bar.

``atmos_top_pressure``
  A ``float`` that contains the top pressure of the atmosphere in bar.

``temperature_profile_model``
  A ``str`` that contains the choice of the temperature profile model. The available options can be found :ref:`here <sec:temperature_profiles>`.

``temperature_profile_parameters``
  A ``list[str]`` that contains the parameters of the temperature profile model.

``radiative_transfer_model``
  A ``str`` that contains the choice of the radiative transfer model. The available options can be found :ref:`here <sec:forward_model_em>`.

``radiative_transfer_parameters``
  A ``list[str]`` that contains the parameters of the radiative transfer model.

``chemistry_model``
  A ``list[str]`` that contains the choice of the chemistry model(s). The available options can be found :ref:`here <sec:chemistry_models>`.

``chemistry_parameters``
  A ``list[list[str]]`` that contains the parameters of the chemistry model(s).

``opacity_species_symbol``
  A ``list[str]`` that contains the symbols/formulas of the opacity species used in the radiative transfer calculations.

``opacity_species_folder``
  A ``list[str]`` that contains the folder where opacity for each species are stored.

The third constructor is similar to the second one but includes additional arguments for cloud models:

.. code:: python

   model_config = bear.EmissionModelConfig(
      int nb_grid_points,
      float atmos_bottom_pressure,
      float atmos_top_pressure,
      str temperature_profile_model,
      list[str] temperature_profile_parameters,
      str radiative_transfer_model,
      list[str] radiative_transfer_parameters,
      list[str] chemistry_model,
      list[list[str]] chemistry_parameters,
      list[str] opacity_species_symbol,
      list[str] opacity_species_folder,
      list[str] cloud_model,
      list[list[str]] cloud_model_parameters)

The additional arguments are:

``cloud_model``
  A ``list[str]`` that contains the choice of the cloud model(s). The available options can be found :ref:`here <sec:cloud_models>`.

``cloud_model_parameters``
  A ``list[list[str]]`` that contains the parameters of the cloud model(s).


PostProcess Config
==================

The configuration class for the emission spectroscopy model post process has the following structure:

.. code:: python

   class EmissionPostProcessConfig:
    __init__(str file_path)
    __init__(
      bool save_temperatures,
      bool save_effective_temperatures,
      bool save_spectra, 
      bool save_contribution_functions,
      list[str] chemical_species_symbols)
    save_temperatures : bool = True
    save_effective_temperatures : bool = True
    save_spectra : bool = True
    save_contribution_functions : bool = False,
    delete_sampler_files : bool = False

The class contains two different constructors. The first one takes a string argument that is the path to the
secondary eclipse model post process configuration file: `file_path`.

.. code:: python

   post_process_config = bear.EmissionPostProcessConfig(file_path)

The second constructor takes the most important parameters of the post process as arguments:

.. code:: python

   post_process_config = bear.EmissionPostProcessConfig(
      bool save_temperatures,
      bool save_effective_temperatures,
      bool save_spectra,
      bool save_contribution_functions,
      list[str] chemical_species_symbols)

The different options are discussed :ref:`here <sec:forward_model_em>`.

.. warning::

   In the current Python wrapper, the *attribute* ``save_effective_temperatures`` of ``bear.EmissionPostProcessConfig``
   is erroneously bound to the ``save_temperatures`` member of the underlying C++ class. Reading or writing
   ``post_process_config.save_effective_temperatures`` therefore accesses ``save_temperatures`` instead. The
   constructor argument is not affected; set this option through the constructor (or the ``post_process.toml`` file)
   and do not modify the attribute afterwards.


Flat Line
*********

The flat-line forward model produces a constant (flat) spectrum and is mainly used as a null model for
model comparison. It does not have a dedicated Python model configuration class; only a post process
configuration class is exposed.

PostProcess Config
==================

The configuration class for the flat-line model post process has the following structure:

.. code:: python

   class FlatLinePostProcessConfig:
    __init__(str file_path)
    __init__(bool save_spectra)
    save_spectra : bool = False
    delete_sampler_files : bool = False

The class contains two different constructors. The first one takes a string argument that is the path to the
flat-line model post process configuration file: `file_path`.

.. code:: python

   post_process_config = bear.FlatLinePostProcessConfig(file_path)

The second constructor takes the most important parameter of the post process as an argument:

.. code:: python

   post_process_config = bear.FlatLinePostProcessConfig(
      bool save_spectra)


Base Class and Direct Forward-Model Bindings
********************************************

All forward-model configuration classes and post process configuration classes described above derive from the
common base class ``bear.GenericConfig``. This class is exposed only so that the constructors of ``bear.Retrieval``
and ``bear.PostProcess`` can accept any of the derived configuration objects; it has no constructor, attributes, or
methods of its own and cannot be instantiated from Python.

Besides the retrieval interface, the wrapper also exposes three forward models directly. They allow the user to
calculate a spectrum for a given pressure-temperature structure and chemical composition without setting up a
retrieval, observations, or priors. The high-level examples of this usage are the scripts
``python/transmission/bear_transmission_spectrum.py``, ``python/secondary_eclipse/bear_secondary_eclipse_spectrum.py``,
and ``python/emission/bear_emission_spectrum.py`` with their Python helper classes in the respective ``src``
subfolders (see :ref:`here <sec:python_interface>`). The raw classes have the following structure:

.. code:: python

   class TransmissionModel:
    __init__(
      bear.Config config,
      bear.SpectralGrid spectral_grid,
      int nb_grid_points,
      list[str] opacity_species_symbol,
      list[str] opacity_species_folder)
    list[float] calcSpectrum(
      float surface_gravity,
      float planet_radius,
      float radius_ratio,
      list[float] pressure,
      list[float] temperature,
      list[str] species_symbol,
      list[list[float]] mixing_ratios,
      list[list[float]] cloud_optical_depth,
      bool use_variable_gravity)

   class OccultationModel:
    __init__(
      bear.Config config,
      bear.SpectralGrid spectral_grid,
      int nb_grid_points,
      list[float] stellar_spectrum_wavelengths,
      list[float] stellar_spectrum_flux,
      list[str] opacity_species_symbol,
      list[str] opacity_species_folder)
    list[float] calcSpectrum(
      float surface_gravity,
      float radius_ratio,
      list[float] pressure,
      list[float] temperature,
      list[str] species_symbol,
      list[list[float]] mixing_ratios,
      list[list[float]] cloud_optical_depth)

   class EmissionModel:
    __init__(
      bear.Config config,
      bear.SpectralGrid spectral_grid,
      int nb_grid_points,
      str radiative_transfer_model,
      list[str] radiative_transfer_parameters,
      list[str] opacity_species_symbol,
      list[str] opacity_species_folder)
    list[float] calcSpectrum(
      float surface_gravity,
      float radius,
      float distance,
      list[float] pressure,
      list[float] temperature,
      list[str] species_symbol,
      list[list[float]] mixing_ratios,
      list[list[float]] cloud_optical_depth,
      bool use_variable_gravity)

The constructors take a ``bear.Config`` object, in which at least ``use_gpu``, ``forward_model_type``,
``cross_section_file_path``, ``spectral_disecretisation``, and ``spectral_resolution`` should be set, and a
``bear.SpectralGrid`` object created with the three-argument constructor, i.e. with an explicit wavelength range. The
argument ``nb_grid_points`` is the number of atmospheric levels; ``opacity_species_symbol`` and
``opacity_species_folder`` describe the opacity sources as in the forward-model configuration classes. The
``OccultationModel`` additionally requires a stellar spectrum given as a list of wavelengths (in micron) and the
corresponding flux values (in W m\ :sup:`-2` micron\ :sup:`-1`), while the ``EmissionModel`` requires the choice of the radiative transfer model and its
parameters (see :ref:`here <sec:forward_model_em>`).

The ``calcSpectrum`` methods compute a single spectrum on the spectral grid and return it as a ``list[float]`` in the
order of the ``wavenumber_list`` of the spectral grid (i.e. in increasing wavenumber, decreasing wavelength). The
arguments ``pressure`` (in bar) and ``temperature`` (in K) are lists with ``nb_grid_points`` entries. The argument
``species_symbol`` is a list of chemical species symbols and ``mixing_ratios`` the corresponding volume mixing ratios,
where the first index refers to the atmospheric level and the second one to the species. The argument
``cloud_optical_depth`` contains the optical depth of each atmospheric *layer* (``nb_grid_points - 1`` entries) at
every wavenumber of the spectral grid. The surface gravity is given in cgs units (cm/s\ :sup:`2`); radii and distances
are given in cm. For the transmission model, ``radius_ratio`` is the ratio of the planet to the stellar radius and
``planet_radius`` the radius at the bottom of the atmosphere. For the emission model, ``radius`` and ``distance`` are
used to scale the emitted flux to the one observed at the given distance. The flag ``use_variable_gravity`` determines
if the gravity is allowed to vary with altitude.
