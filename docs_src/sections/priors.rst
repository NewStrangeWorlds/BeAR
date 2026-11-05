.. _sec:prior_distributions:

Prior distributions
===================

Each free parameter of a retrieval calculations needs to have an associated prior distribution. The prior distribution
describes any prior knowledge that is available about the parameter. The number and type of the free parameters
depend on the chosen forward model and user-specified configuration options.

BeAR supports the following types of distributions:

  - ``delta`` - a delta distribution with a fixed value

  - ``uniform`` - a uniform distribution with a lower and upper bound

  - ``log_uniform`` - a log-uniform distribution with a lower and upper bound

  - ``gaussian`` - a Gaussian distribution with a mean and standard deviation

  - ``linked`` - links this prior to that of another parameter

A special case is the ``linked`` distribution. This distribution links the prior distribution
of one parameter to that of another. Thus, these two parameters will always have the same value during a retrieval
calculation.


Prior distributions units
.........................

BeAR also supports units for its prior distributions. The following units are currently taken into account:

  - ``Rs`` or ``Rsun`` - the solar radius

  - ``Rj`` or ``Rjupiter`` - Jupiter's radius

  - ``Re`` or ``Rearth`` - Earth's radius

  - ``Ms`` or ``Msun`` - the solar mass

  - ``Mj`` or ``Mjupiter`` - Jupiter's mass

  - ``Me`` or ``Mearth`` - Earth's mass

  - ``pc`` - distance in parsec

  - ``ly`` - distance in light years

Priors without units are assumed to be in cgs units.


Prior distributions file
........................

The ``priors.config`` file contains the information on the prior distributions of the free parameters.
The file has the following structure:

.. include:: ../examples/priors_example1.config
   :literal:

Each row describes the prior for one free parameter. The first column lists the distribution type of the prior,
the second column the model parameter name, and the remaining columns the parameters of the distribution, while
the optional, last column is the unit of the parameter. The type and number of distribution parameters depend on
the chosen distribution type.

The rows are matched to the free parameters of the forward model **by name**: the parameter name in the second
column must exactly match one of the model's parameter names. The order of the rows in the file therefore does
not matter. Every model parameter must have exactly one matching row. Missing parameter names, names that do not
correspond to any model parameter, and duplicate names are all treated as hard errors and will stop the retrieval.

A special case is the ``linked`` distribution. This distribution links the prior distribution
of one parameter to that of another, so that the two parameters always share the same value during a retrieval.
The row has the form ``linked <this_param_name> <target_param_name>``, i.e. the second column is the parameter
being defined and the third column is the **name** of the target parameter it should be linked to.
An example of this is shown below:

.. include:: ../examples/priors_example2.config
   :literal:

Here, the prior for the CO2 mixing ratio (``chem_co2_mr``) is linked to that of the methane mixing ratio
(``chem_ch4_mr``). Thus, CO2 will always have the same mixing ratio as methane for this retrieval setup. A prior
may not link to itself and may not link to another parameter that is itself a ``linked`` prior; both cases are
rejected. It is also important to note that BeAR cannot check the consistency of the linked parameters. For
example, if the linked target were a temperature, the resulting mixing ratios of CO2 would make no sense. It is
the user's responsibility to ensure that the linked parameters are consistent.


Required parameter names
........................

Because the priors are matched by name, the ``priors.config`` file must contain a row for each free parameter of
the chosen forward model, in any order. The set of names that must be supplied is, grouped by their origin:

  - general forward model parameters, as discussed in the :ref:`section <sec:forward_models>` on forward models

  - parameters for the chosen chemistry models, as discussed in the :ref:`section <sec:chemistry_models>` on chemistry models

  - temperature profile parameters, as discussed in this :ref:`section <sec:temperature_profiles>`

  - cloud model parameters, as discussed in the :ref:`section <sec:cloud_models>` on cloud models

  - parameters for optional modules that can be used for a specific forward model

  - observational offset parameters, as discussed in the :ref:`section <sec:observational_data>` on observations

A general prior configuration file therefore contains the following parameters (shown here in the canonical
grouping, but the rows may be listed in any order):

.. include:: ../examples/priors_example3.config
   :literal:

Two sets of parameters are handled automatically and their prior behaviour depends on the retrieval configuration:

  - When ``use_error_inflation`` is set to ``true`` in ``retrieval.toml``, BeAR automatically appends an internal
    error-inflation prior (a uniform ``error exponent`` prior whose bounds are derived from the observational
    uncertainties). This prior must **not** be added to ``priors.config`` by the user.

  - When high-resolution observations are present, the retrieval appends a high-resolution tail of parameters named
    ``kp``, ``vsys``, ``dphi`` and, optionally, ``alpha``. Priors for ``kp``, ``vsys`` and ``dphi`` must be supplied
    in ``priors.config``; a prior for ``alpha`` is only added as a free parameter if a row named ``alpha`` is present.
    Full documentation of the high-resolution retrieval will be added later.
