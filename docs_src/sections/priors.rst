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

  - ``pc`` - distance in parsec

  - ``ly`` - distance in light years

Priors without units are assumed to be in cgs units.


Prior distributions file
........................

The ``priors.config`` file contains the information on the prior distributions of the free parameters.
The file has the following structure:

.. include:: ../examples/priors_example1.config
   :literal:

The first column lists the distribution type of the prior, the second column the model parameter name, and the 
remaining columns the parameters of the distribution, while the optional, last column is the unit of the parameter.

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


Order of prior distributions
............................

The order of the prior distributions in the ``priors.config`` file is very important. In general they have to appear in the 
following order:

  - general forward model parameters, as discussed in the :ref:`section <sec:forward_models>` on forward models

  - priors for the chosen chemistry models, as discussed in the :ref:`section <sec:chemistry_models>` on chemistry models

  - temperature profile priors, as discussed in this :ref:`section <sec:temperature_profiles>`

  - cloud model priors as discussed in the :ref:`section <sec:cloud_models>` on cloud models

  - priors for optional modules that can be used for a specific forward model

  - observational offset priors as discussed in the :ref:`section <sec:observational_data>` on observations

A general prior configuration file, therefore has the following structure:

.. include:: ../examples/priors_example3.config
   :literal:
