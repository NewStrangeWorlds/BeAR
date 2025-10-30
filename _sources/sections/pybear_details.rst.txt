
.. _sec:pybear_details:

##############################
Detailed Description of pyBeAR
##############################


Retrieval and PostProcess classes
#################################

These two classes are the main classes of the pyBeAR package. The Retrieval class is used to set up and perform
the retrieval calculations, while the PostProcess class is used to perform post retrieval calculations based
on the posterior samples obtained from a retrieval run.


pybear.Retrieval class
**********************

This class describes the retrieval object that is used to perform the retrieval calculations. It is mainly responsible for converting
the model parameters into physical values, computing the forward model spectra, and calculating the likelihoods. likelihood calculations can
either done by the retrieval itself using the internal MultiNest sampler and likelihood functions of BeAR or by external functions written in Python.
For the latter, this class can return the model spectra for a given set of parameters that can then be compared to the observational data in Python.

The class has the following structure:

.. code:: python

   class Retrieval
    __init__(
      pyBear.Config config)
    __init__(
      pyBear.Config config, 
      pyBear.ModelConfig model_config, 
      list[pyBear.Observation] observations,
      list[pyBear.Prior] priors)
    pyBear.SpectralGrid spectral_grid
    run()
    int nbParameters()
    pair(list[float], list[float]) convertCubeParameters(
      list[float] cube_parameters)
    list[float] convertToPhysicalParameters(
      list[float] cube_parameters)
    float computeLikelihood(
      list[float] physical_parameters)
    pyBear.ForwardModelOutput computeModel(
      list[float] physical_parameters,
      bool return_high_res_spectrum)
    pyBear.AtmosphereOutput computeAtmosphereStructure(
        list[float] physical_parameters,
        list[str] species_symbols)

The class contains two different constructors. The first one takes a single argument that is the general configuration
of BeAR:

.. code:: python

   retrieval = pybear.Retrieval(config)

In this case, the other configuration files for the chosen forward model, the observational data, and the prior information
are read in during the setup of the object from the folder set in the ``config.retrieval_folder_path`` variable.

The second constructor takes additional arguments for the model configuration, observations, and priors:

.. code:: python

   retrieval = pybear.Retrieval(
    config,
    model_config,
    observations,
    priors)

These config objects need to be created beforehand with the corresponding classes of pyBeAR as described in the following sections.

The Retrieval class has the following attributes:

``spectral_grid``
  A ``pyBear.SpectralGrid`` object that contains the spectral grid information for the retrieval.

The Retrieval class has the following methods:

``run()``
  This method starts the retrieval calculation using the internal MultiNest sampler and likelihood functions of BeAR.

``nbParameters()``
  This method returns an ``int`` that contains the number of free parameters of the retrieval.

``convertCubeParameters(cube_parameters)``
  This method takes a ``list[float]`` of hypercube parameters as input and returns a ``pair(list[float], list[float])``
  that contains the values of the hypercube parameters into their prior distributions and their physical values. 
  For the latter, the units of the priors are taken into account.

``convertToPhysicalParameters(cube_parameters)``
  This method takes a ``list[float]`` of prior parameters as input and returns a ``list[float]`` that contains
  the corresponding physical parameters.

``computeLikelihood(physical_parameters)``
  This method takes a ``list[float]`` of physical parameters as input, internally calculates a forward model for these parameters, 
  compares it to the observational data and returns a ``float`` that contains the corresponding likelihood value for this parameter set.

``computeModel(physical_parameters, return_high_res_spectrum)``
  This method takes a ``list[float]`` of physical parameters as input and returns an object of type ``pyBear.ForwardModelOutput`` that contains
  the corresponding model spectra. If the argument ``return_high_res_spectrum`` is set to ``True``, the high-resolution spectrum
  will also be returned in the output object.

``computeAtmosphereStructure(physical_parameters, species_symbols)``
  This method takes a ``list[float]`` of physical parameters as input and returns a ``pyBear.AtmosphereOutput`` object that contains
  the corresponding atmosphere structure for these parameters. The argument ``species_symbols`` is a ``list[str]`` that contains the chemical 
  symbols of the species that should be included in the atmosphere structure output.


pybear.PostProcess class
************************

The postprocess class is used to perform post retrieval calculations based on the posterior samples obtained from a retrieval run.
This includes the calculation of model spectra for the posterior samples, effective temperatures, and other derived parameters. 
This class will either read in the posterior sample or take them as an argument. The class is a child class of the ``pyBear.Retrieval`` 
class and thus inherits most of its attributes and methods.

The class has the following structure, where only the additional constructors and methods are listed:

.. code:: python

    class PostProcess(pyBear.Retrieval):
      __init__(
        pyBear.Config config)
      __init__(
        pyBear.Config config, 
        pyBear.ModelConfig model_config, 
        pyBear.PostProcessConfig post_process_config,
        list[pyBear.ObservationInput] observations,
        list[pyBear.Prior] priors)
      run()
      run(str posterior_file_path)
      run(
        list[list[float]] &posteriors,
        list[float] &log_likelihoods)

The constructors are similar to those of the ``pyBear.Retrieval`` class except that the second one takes an
additional argument for the postprocess configuration. 

The PostProcess class has the following additional methods:

``run()``
  This method starts the post process calculations using the standard MultiNest output files located in the folder 
  set in the ``config.retrieval_folder_path`` variable. 
  The output will be saved in the folder set by the ``config.post_output_path`` variable.

``run(output_folder)``
  This method takes a ``str`` argument that is the path to the output posterior sample file. 
  The output will be saved in the folder set by the ``config.post_output_path`` variable.

``run(posteriors, log_likelihoods)``
  This method takes a ``list[list[float]]`` of posterior samples and a ``list[float]`` of corresponding log-likelihoods
  as arguments and starts the post process calculations based on these samples.
  The output will be saved in the folder set by the ``config.post_output_path`` variable.


pybear.Config
#############

This class describes the general model configuration of BeAR and is equivalent to the information located in a 
``retrieval.config`` file. The class is used to load the configuration files of BeAR and to set up the retrieval object.

It has the following attributes and methods:

.. code:: python

   class Config:
       __init__()
       __init__(str file_path)
       __init__(
          bool use_gpu, 
          str forward_model_type, 
          str opacity_path,
          str spectral_discretisation,
          float resolution,
          str multinest_output_folder,
          str post_process_output_folder)
      forward_model_type : str = ""
      retrieval_folder_path : str = ""
      multinest_output_path : str = ""
      post_output_path : str = ""
      wavenumber_file_path : str = ""
      cross_section_file_path : str = ""
      spectral_disecretisation : str = ""
      spectral_resolution : float = 0
      multinest_print_iter_values : bool = False
      multinest_efficiency : float = 0.8
      multinest_nb_living_points : int = 800
      multinest_nb_iterations : int = 0
      multinest_feedback : bool = False
      use_error_inflation : bool = False
      use_gpu : bool = False
      nb_omp_processes : int = 0

The class contains three different constructors. The first one is the default constructor that does not take any arguments.

.. code:: python

   model_config = pybear.Config()

In this case the object will be created will all values set to their default values. The second constructor takes a string
argument that is the path to the retrieval folder. This is the folder where the configuration file ``retrieval.config`` is 
located.

.. code:: python

   model_config = pybear.Config(folder)

The function will load the configuration file and set the attributes of the object accordingly. The third constructor is the
most flexible one and allows the user to set the most important attributes of the object directly. The arguments are:

.. code:: python

   model_config = pybear.Config(
    bool use_gpu, 
    str forward_model_type, 
    str cross_section_file_path,
    str spectral_discretisation,
    float spectral_resolution,
    str multinest_output_folder,
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

``multinest_output_path``
  A ``str`` that contains the path to the folder where the MultiNest output files will be stored.

``post_output_path``
  A ``str`` that contains the path to the folder where the post process output files will be stored.

``wavenumber_file_path``
  A ``str`` that contains the path to the file that contains the wavenumber grid. If the standard HELIOS-k opacities 
  with a constant wavenumber step of 0.01 cm^-1 are used, this file is not needed.

``cross_section_file_path``
  A ``str`` that contains the path to the file that contains the opacity data.

``spectral_discretisation``
  A ``str`` that contains the choice of the spectral discretisation. The available options can be found :ref:`here <sec:config_files>`

``spectral_resolution``
  A ``float`` that contains the spectral resolution or constant wavelength/wavenumber step of the data.

``multinest_print_iter_values``
  A ``bool`` that determines if the model parameters and iteration values of MultiNest are printed to the terminal.

``multinest_efficiency``
  A ``float`` that contains the efficiency of the MultiNest sampler.

``multinest_nb_living_points``
  An ``int`` that contains the number of living points of the MultiNest sampler.

``multinest_nb_iterations``
  An ``int`` that contains the number of iterations of the MultiNest sampler. A value of 0 implies that no fixed number of iterations is 
  set and MultiNest will iterate until it converges.

``multinest_feedback``
  A ``bool`` that determines if the feedback of the MultiNest sampler is enabled. This includes printing regular 
  updates of the iteration to the terminal.

``use_error_inflation``
  A ``bool`` that determines if the error inflation is enabled.

``use_gpu``
  A ``bool`` that determines if the GPU is used for the calculations.

``nb_omp_processes``
  An ``int`` that contains the number of OpenMP processes that are used for the calculations. 
  A number of 0 implies that the maximum available number of processes will be used.


pybear.SpectralGrid
###################

This class describes the spectral grid that is used for the retrieval calculations. The class has the following structure:

.. code:: python
  
   class SpectralGrid
    __init__(
      pyBear.GlobalConfig config, 
      float min_wavenumber, 
      float max_wavenumber)
    wavenumber_list : list[float]
    wavelength_list : list[float]
    float wavelengthToWavenumber(
      float wavelength)

The constructor takes the general config as an argument as well as the minimum and maximum wavenumber of the spectral grid:

.. code:: python

   spectral_grid = pybear.SpectralGrid(
    config,
    min_wavenumber,
    max_wavenumber)

The spectral grid will be created based on the provided wavenumber range and the spectral discretisation and resolution
set in the configuration. Note that the retrieval and post process classes will automatically create the spectral grid during their setup,
so the user usually does not need to create this object directly.

The SpectralGrid class has the following attributes and methods:

``wavenumber_list``
  A ``list[float]`` that contains the wavenumber grid in units of cm^-1.

``wavelength_list``
  A ``list[float]`` that contains the wavelength grid in units of micron.

``float wavelengthToWavenumber(float wavelength)``
  This method takes a ``float`` argument that is a wavelength in units of micron and returns a ``float``
  that contains the corresponding wavenumber in units of cm^-1.


pybear.Observation
##################

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
      std type, 
      list[float] wavelengths,
      list[float] observation_data,
      list[float] observation_error)
    __init__(
      str name, 
      std type, 
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

   observation = pybear.Observation()

In this case the object will be created will all values set to their default values. The second constructor takes two string
argumentse that is the name of the observation and its type:

.. code:: python

   observation = pybear.Observation(name, type)

The third constructor is used for spectroscopy observations that are given as a list of wavelengths, observation data, 
and observation errors:

.. code:: python

   observation = pybear.Observation(
    name, 
    type, 
    wavelengths, observation_data, 
    observation_error)

Finally, the fourth constructor is used for photometry or band-spectroscopy observations that are given as a list of wavelength bins, observation data, and observation errors:

.. code:: python

   observation = pybear.Observation(
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


pybear.Prior
############

This class describes a specific model prior distribution used in the retrieval as described :ref:`here <sec:prior_distributions>`.
The class has the following structure:

.. code:: python

   class Prior
    __init__(
      str description, 
      str type, 
      list[float] parameters)
    __init__(
      str description, 
      str type, 
      list[float] parameters,
      str unit)
    description : str = ""
    type : str = ""
    parameter : list[float]
    unit : str = ""

The class contains two different constructors. The first one is used for priors that do not have a unit:

.. code:: python

   prior = pybear.Prior(
    description, 
    type, 
    parameters)

For example, a uniform prior can be created via:

.. code:: python

   prior = pybear.Prior(
    "log g", 
    "uniform", 
    [2.0, 5.0]) 

The second constructor is used for priors that have a unit:

.. code:: python

   prior = pybear.Prior(
    description, 
    type, 
    parameters,
    unit)

For example, a Gaussian prior in solar radii can be created via:

.. code:: python

   prior = pybear.Prior(
    "stellar_radius", 
    "gaussian", 
    [1.5, 0.1],
    "Rs")

All attributes can also be set directly after the object has been created. For example, the unit can be changed via:

.. code:: python

   prior.unit = "Rj"

The Prior class has the following attributes:

``description``
  A ``str`` that contains the description/name of the prior.

``type``
  A ``str`` that contains the type of the prior. The available options can be found :ref:`here <sec:prior_distributions>`

``parameters``
  A ``list[float]`` that contains the parameters of the prior. For a prior that only requires one parameter, 
  this variable still needs to be a list.

``unit``
  A ``str`` that contains the optional unit of the prior. A list of possible units can be found :ref:`here <sec:prior_distributions>`





pybear.ForwardModelOutput
#########################

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
  via the ``pyBear.Observation`` objects.

``spectrum``
  A ``list[float]`` that contains the high-resolution spectrum calculated on the spectral grid of the retrieval.


pybear.AtmosphereOutput
#######################

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


*******************************************
Forward Model Config and PostProcess Config
*******************************************

Every forward model in BeAR has its own configuration that is used to set up the model and an additional one 
for the post process calculations. The configuration parameters are different for every model type and are 
described in the :ref:`section <sec:forward_models>` on forward models. Below, the Python classes for these
configurations are described.

Transmission Spectrum
#####################


Model Config
************

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
      list[str] opacity_species_folder
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
    cloud_model_parameters : list[str]
    modules : list[str]
    modules_parameters : list[list[str]]
    opacity_species_symbol : list[str]
    opacity_species_folder : list[str]

The class contains three different constructors. The first one takes a string argument that is the path to the
transmission model configuration file: `file_path`. 

.. code:: python

   model_config = pybear.TransmissionModelConfig(file_path)

The second constructor takes the most important parameters of the model as arguments:

.. code:: python

   model_config = pybear.TransmissionModelConfig(
      int nb_grid_points,
      float atmos_bottom_pressure,
      float atmos_top_pressure,
      str temperature_profile_model,
      list[str] temperature_profile_parameters,
      list[str] chemistry_model,
      list[list[str]] chemistry_parameters,
      list[str] cloud_model,
      list[list[str]] cloud_model_parameters)

The constructor has the following arguments:

``nb_grid_points``
  An ``int`` that contains the number of grid points of the atmosphere.

``atmos_bottom_pressure_``
  A ``float`` that contains the bottom pressure of the atmosphere in bar.

``atmos_top_pressure_``
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

``opacity_species_folder``
  A ``list[str]`` that contains the folder where opacity for each species are stored.

The third constructor is the most flexible one and takes additional arguments for advanced model options:

.. code:: python

   model_config = pybear.TransmissionModelConfig(
      bool fit_mean_molecular_weight,
      bool fit_scale_height,
      bool use_variable_gravity,
      int nb_grid_points,
      float atmos_bottom_pressure_,
      float atmos_top_pressure_,
      str temperature_profile_model,
      list[str] temperature_profile_parameters,
      list[str] chemistry_model,
      list[list[str]] chemistry_parameters,
      list[str] opacity_species_symbol,
      list[str] opacity_species_folder
      list[str] cloud_model,
      list[str] cloud_model_parameters,
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
******************

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

   post_process_config = pybear.TransmissionPostProcessConfig(file_path)

The second constructor takes the most important parameters of the post process as arguments:

.. code:: python

   post_process_config = pybear.TransmissionPostProcessConfig(
      bool save_temperatures,
      bool save_spectra,
      list[str] chemical_species_symbols)

The different options are discussed :ref:`here <sec:forward_model_transmission>`.


Secondary Eclipse Spectrum
##########################

Model Config
************

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
      list[str] opacity_species_folder
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

   model_config = pybear.OccultationConfig(file_path)

The second constructor takes the most important parameters of the model as arguments:

.. code:: python

   model_config = pybear.OccultationConfig(
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
      list[str] opacity_species_folder
      str stellar_spectrum_model,
      list[str] stellar_model_parameters,
      list[str] cloud_model,
      list[list[str]] cloud_model_parameters)

The constructor has the following arguments:

``nb_grid_points``
  An ``int`` that contains the number of grid points of the atmosphere.

``atmos_bottom_pressure_``
  A ``float`` that contains the bottom pressure of the atmosphere in bar.

``atmos_top_pressure_``
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
  A ``list[str]`` that contains the symbols/formulas of the opacity species used in the radiative transfer calculations.§

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
******************

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

   post_process_config = pybear.OccultationPostProcessConfig(file_path)

The second constructor takes the most important parameters of the post process as arguments:

.. code:: python

   post_process_config = pybear.OccultationPostProcessConfig(
      bool save_temperatures,
      bool save_spectra,
      bool save_contribution_functions,
      list[str] chemical_species_symbols)

The different options are discussed :ref:`here <sec:forward_model_se>`.


Secondary Eclipse Spectrum with Planetary Blackbody
###################################################

Model Config
************

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

   model_config = pybear.OccultationBlackBodyConfig(file_path)

The second constructor takes the most important parameters of the model as arguments:

.. code:: python

   model_config = pybear.OccultationBlackBodyConfig(
      str stellar_spectrum_model,
      list[str] stellar_model_parameters)

The constructor has the following arguments:

``stellar_spectrum_model``
  A ``str`` that contains the choice of the stellar spectrum model. The available options can be found :ref:`here <sec:stellar_spectra>`.

``stellar_model_parameters``
  A ``list[str]`` that contains the parameters of the stellar spectrum model.

PostProcess Config
******************

The configuration class for the secondary eclipse model post process has the following structure:

.. code:: python

   class OccultationBlackBodyPostProcessConfig:
    __init__(str file_path)
    __init__(bool save_spectra,)
    save_spectra : bool = True
    delete_sampler_files : bool = False

The class contains two different constructors. The first one takes a string argument that is the path to the
secondary eclipse model post process configuration file: `file_path`.

.. code:: python

   post_process_config = pybear.OccultationPostProcessConfig(file_path)

The second constructor takes the most important parameters of the post process as arguments:

.. code:: python

   post_process_config = pybear.OccultationPostProcessConfig(
      bool save_spectra)

The different options are discussed :ref:`here <sec:forward_model_se_bb>`.


Emission Spectrum
#################

Model Config
************

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
      str chemistry_model,
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
      list[str] opacity_species_folder
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

   model_config = pybear.EmissionModelConfig(file_path)

The second constructor takes the most important parameters of the model as arguments:

.. code:: python

   model_config = pybear.EmissionModelConfig(
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

``atmos_bottom_pressure_``
  A ``float`` that contains the bottom pressure of the atmosphere in bar.

``atmos_top_pressure_``
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
  A ``list[str]`` that contains the symbols/formulas of the opacity species used in the radiative transfer calculations.§

``opacity_species_folder``
  A ``list[str]`` that contains the folder where opacity for each species are stored.

The third constructor is similar to the second one but includes additional arguments for cloud models:

.. code:: python

   model_config = pybear.EmissionModelConfig(
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
      list[str] opacity_species_folder
      list[str] cloud_model,
      list[list[str]] cloud_model_parameters)

The additional arguments are:

``cloud_model``
  A ``list[str]`` that contains the choice of the cloud model(s). The available options can be found :ref:`here <sec:cloud_models>`.

``cloud_model_parameters``
  A ``list[list[str]]`` that contains the parameters of the cloud model(s).


PostProcess Config
******************

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

   post_process_config = pybear.EmissionPostProcessConfig(file_path)

The second constructor takes the most important parameters of the post process as arguments:

.. code:: python

   post_process_config = pybear.EmissionPostProcessConfig(
      bool save_temperatures,
      bool save_effective_temperatures,
      bool save_spectra,
      bool save_contribution_functions,
      list[str] chemical_species_symbols)

The different options are discussed :ref:`here <sec:forward_model_em>`.