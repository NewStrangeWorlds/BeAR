
.. _sec:forward_models:

Forward Models
==============

BeAR currently includes the following forward models:

  - :ref:`Transmission spectrum <sec:forward_model_transmission>`

  - :ref:`Emission spectrum <sec:forward_model_em>`

  - :ref:`Secondary eclipse spectrum <sec:forward_model_se>`

  - :ref:`Secondary eclipse spectra with a planetary blackbody <sec:forward_model_se_bb>`

  - :ref:`Flat line <sec:forward_model_flat>`

The latter two models are usually used to check if the observational
data is more likely explained by more simpler models, such as a flat line.

Each model requires its own config file ``forward_model.config`` that are 
discussed below. 


Model postprocessing
--------------------

After the retrieval calculations are finished, the  BeAR will perform a postprocessing on 
the resulting posterior sample. Depending on the chosen forward model, different 
postprocessing steps can be used, including saving all posterior spectra or temperature 
profiles. Details on the output files can be found in the :ref:`section <sec:output_files>` 
on output files.

These postprocessing steps can also be configured in the optional
``post_process.config`` file. If this file is not present, BeAR will use 
default settings for the postprocessing. Depending on the forward model, 
available options are currently:

  - ``Delete unused MultiNest files`` - This will option will delete all MultiNest files that
    are not used in the post processing. Available options: ``Yes`` or ``No``.

  - ``Save spectra`` - Compute the spectra for each model in the posterior sample. Spectral
    will be saved for each observational/instrument individually. Additionally, a high-resolution
    spectrum of the best-fit model will be computed and saved. Available options: ``Yes`` or ``No``.

  - ``Save temperature structures`` - Compute and save the temperature profile for each model in the posterior
    sample. Available options: ``Yes`` or ``No``.

  - ``Save effective temperatures`` - Compute and save the effective temperature for each model in posterior sample.
    Available options: ``Yes`` or ``No``.

  - ``Save contribution functions`` - Compute and save the contribution functions for each observational/instrument 
    for the best-fit model. Available options: ``Yes`` or ``No``.

  - ``Save chemical species`` - Save the mixing ratio profiles for selected chemical species for all 
    posterior samples. The options for this paramter are the forumulas of the chemical species that should be saved,
    separated by white spaces. For example, ``H2O CO2`` will save the mixing ratios of water and carbon dioxide.
    Species that BeAR does not know will be ignored. 



.. _sec:forward_model_flat:

Flat line
---------

This model simply fits a flat line through the observational data. Usually, this
model is used to test if the use of a more complex model is warranted to explain
the observation. The test is normally done by comparing their 
Bayesian evidences.

In the main retrieval config file ``retrieval.config`` it is selected by choosing:

.. code:: 

   #Forward model
   flat_line


Model config file
.................

The flat line model does not need a ``forward_model.config`` file since it
has no configurable parameters. In the prior distribution file, this model 
requires a single free parameter that determines the flat line. 

This parameter needs to have the same units as the observational data. For example,
in case of a transmission spectrum, this parameter refers to the transit depth in ppm, while
for an emission spectrum, it needs to have units of 
:math:`\mathrm{W} \mathrm{m^{-2}} \mathrm{\mu m^{-1}}`.


Model postprocessing
....................

The optional postprocessing file ``post_process.config`` has the following structure and
options:

.. include:: ../examples/post_process_flatline.config
   :literal:

The default options that are used when this file is not present are:

  - ``Delete unused MultiNest files`` : No

  - ``Compute spectra`` : No



.. _sec:forward_model_transmission:

Transmission spectrum
----------------------

This forward model computes the wavelength-dependent transit depth 
:math:`D(\lambda)` of an exoplanet atmosphere, given by

.. math::
  D(\lambda) = \left(\frac{R_p(\lambda)}{R_*}\right)^2 \ ,

where :math:`R_p(\lambda)` is the wavelength-dependent planetary radius
and :math:`R_*` the radius of the host star. In BeAR, :math:`D(\lambda)`
has units of ppm.

In the retrieval config file ``retrieval.config`` it is selected by choosing:

.. code:: 

   #Forward model
   transmission

The transmission spectrum forward model has three general free parameters that
have to be added to the priors configuration file in the following order:

  - logarithm of surface gravity :math:`\log g` in cgs units

  - planet radius

  - stellar radius


Model config file
.................

The ``forward_model.config`` file for the transmission spectrum model has
the following structure:

.. include:: ../examples/forward_model_transmission.dat
   :literal:

The first three entries refer to the vertical discretisation of the atmosphere.
This includes the number of atmospheric layers as well as the bottom and 
top-of-atmosphere pressures in units of bars.

With the next entry, the retrieval of the mean molecular or the scale height
can be chosen. Usually, BeAR will determine the mean molecular based on
the chemical abundances of the species included in the retrieval.
Sometimes, however, the background species in an atmosphere might not be
known. For such a case, BeAR can use the mean molecular weight as a free 
parameter. 

Furthermore, for transmission spectra, the surface gravity, the mean molecular
weight, and the temperature might become degenerate in a retrieval. This is
often the case when no constraints on the surface gravity can be provided
and the dominating background specise is not known. In such a scenario, BeAR
can use the (constant) atmospheric scale height in units of km as a free 
parameter.

For the standard case, the keyword ``no`` needs to be used here. If BeAR should
use the mean molecular weight as a free parameter ``mmw`` is used, while
``sh`` is used when the scale should be used as a free parameter instead. If either
the mean molecular weight or the scale height is chosen as a free parameter,
a corresponding prior needs to be added as a fourth model parameter in the
prior configuration file.

With the next two entries, the parametrisation for the 
:ref:`temperature profile <sec:temperature_profiles>` and
the :ref:`cloud models <sec:cloud_models>` are set. After that, optional
modules can be added to the forward model.

Finally, the different :ref:`chemistry models <sec:chemistry_models>`, chemical species, and the
:ref:`opacity sources <sec:opacity_data>` are selected. It is important to note that
chemical species that are used as part of the chemistry models, should normally also have
an associated opacity source. Otherwise, the impact of that species on the 
resulting spectrum might be negligible and, therefore, its abundance will be 
likely be unconstrained.

On the other hand, it is theoretically also possible to add opacity species
without including this species in any of the chemistry models. In this case, the
abundance of this species will be zero and, thus, won't show up in the spectrum.


Model postprocessing
....................

The optional postprocessing file ``post_process.config`` has the following structure and
options:

.. include:: ../examples/post_process_transmission.config
   :literal:

The default options that are used when this file is not present are:

  - ``Delete unused MultiNest files`` : No

  - ``Compute spectra`` : Yes

  - ``Save temperature structures`` : No

  - ``Save chemical species`` : None



.. _sec:forward_model_se:

Secondary eclipse spectrum
---------------------------

This forward model computes the wavelength-dependent secondary eclipse (or occulation) 
depth :math:`D(\lambda)` of an exoplanet atmosphere, given by

.. math::
  D(\lambda) = \frac{F_p(\lambda)}{F_*(\lambda)} \left(\frac{R_p}{R_*}\right)^2 \ ,

where :math:`F_p(\lambda)` is the outgoing flux at top of the planet's atmosphere,
:math:`F_*(\lambda)` is the stellar photospheric flux, :math:`R_p` 
is the planetary radius and :math:`R_*` the radius of the host star. 
In BeAR, :math:`D(\lambda)` has units of ppm.

In the retrieval config file ``retrieval.config`` it is selected by choosing:

.. code:: 

   #Forward model
   secondary_eclipse

The secondary eclipse spectrum forward model has two general free parameters that
have to be added to the priors configuration file in the following order:

  - logarithm of surface gravity :math:`\log g` in cgs units

  - ratio of the planet's and stellar radius :math:`\mathrm{R_p/R_*}`


Model config file
.................

The ``forward_model.config`` file for the secondary eclipse spectrum model has
the following structure:

.. include:: ../examples/forward_model_se.dat
   :literal:

The first three entries refer to the vertical discretisation of the atmosphere.
This includes the number of atmospheric layers as well as the bottom and 
top-of-atmosphere pressures in units of bars.

With the next two entries, the parametrisation for the :ref:`temperature profile <sec:temperature_profiles>`, 
the :ref:`stellar spectrum <sec:stellar_spectra>`, and the :ref:`cloud models <sec:cloud_models>` are set. 

After that, the radiative transfer scheme is chosen. There are currently two 
different options:

  - ``scm`` - the short characteristic method, available for CPU and GPU

  - ``disort`` - the discrete-ordinate solver DISORT, only available for CPU

Finally, the different :ref:`chemistry models <sec:chemistry_models>`, chemical species, and the
:ref:`opacity sources <sec:opacity_data>` are selected. It is important to note that
chemical species that are used as part of the chemistry models, should normally also have
an associated opacity source. Otherwise, the impact of that species on the 
resulting spectrum might be negligible and, therefore, its abundance will be 
likely be unconstrained.

On the other hand, it is theoretically also possible to add opacity species
without including this species in any of the chemistry models. In this case, the
abundance of this species will be zero and, thus, won't show up in the spectrum.


Model postprocessing
....................

The optional postprocessing file ``post_process.config`` has the following structure and
options:

.. include:: ../examples/post_process_se.config
   :literal:

The default options that are used when this file is not present are:

  - ``Delete unused MultiNest files`` : No

  - ``Compute spectra`` : Yes

  - ``Save temperature structures`` : Yes

  - ``Save contribution functions`` : No

  - ``Save chemical species`` : None



.. _sec:forward_model_se_bb:

Secondary eclipse spectrum with planetary blackbody
---------------------------------------------------

This forward model computes the wavelength-dependent secondary eclipse (or occulation) 
depth :math:`D(\lambda)` of an exoplanet atmosphere, given by

.. math::
  D(\lambda) = \frac{F_p(\lambda)}{F_*(\lambda)} \left(\frac{R_p}{R_*}\right)^2 \ ,

where :math:`F_p(\lambda)` is the outgoing flux at top of the planet's atmosphere,
:math:`F_*(\lambda)` is the stellar photospheric flux, :math:`R_p` 
is the planetary radius and :math:`R_*` the radius of the host star. 
In BeAR, :math:`D(\lambda)` has units of ppm.

This model is a special case of the secondary eclipse spectrum model, where the planet's flux
is assumed to be blackbody radiation. This is often used to test if the observational data warrants 
a more complex model or when only a few, usually photometric, data points are available.

In the retrieval config file ``retrieval.config`` it is selected by choosing:

.. code:: 

   #Forward model
   secondary_eclipse_bb

This forward model has two general free parameters that
have to be added to the priors configuration file in the following order:

  - the planet's effective temperature in Kelvin

  - ratio of the planet's and stellar radius :math:`\mathrm{R_p/R_*}`


Model config file
.................

The ``forward_model.config`` file for the secondary eclipse spectrum model has
the following structure:

.. include:: ../examples/forward_model_se_bb.dat
   :literal:

It contains a single option related to description of the stellar spectrum. Information
on the available options can be found in the :ref:`section <sec:stellar_spectra>` on stellar spectra.


Model postprocessing
....................

The optional postprocessing file ``post_process.config`` has the following structure and
options:

.. include:: ../examples/post_process_sebb.config
   :literal:

The default options that are used when this file is not present are:

  - ``Delete unused MultiNest files`` : No

  - ``Compute spectra`` : Yes



.. _sec:forward_model_em:

Emission spectrum
-----------------

This forward model computes the emission spectrum :math:`F(\lambda)` of an exoplanet or brown dwarf atmosphere.
In BeAR, :math:`F(\lambda)` has units of :math:`\mathrm{W} \mathrm{m^{-2}} \mathrm{\mu m^{-1}}`.

In the retrieval config file ``retrieval.config`` it is selected by choosing:

.. code:: 

   #Forward model
   emission

The emission spectrum forward model has two general free parameters that
have to be added to the priors configuration file in the following order:

  - logarithm of surface gravity :math:`\log g` in cgs units

  - a scaling factor :math:`f` for the radius/distance relationship, 
    where the radius is internally set to 1 Jupiter radius

  - distance to the object


Model config file
.................

The ``forward_model.config`` file for the emission spectrum model has
the following structure:

.. include:: ../examples/forward_model_em.dat
   :literal:

The first three entries refer to the vertical discretisation of the atmosphere.
This includes the number of atmospheric layers as well as the bottom and 
top-of-atmosphere pressures in units of bars.

With the next two entries, the parametrisation for the 
:ref:`temperature profile <sec:temperature_profiles>` and
the :ref:`cloud models <sec:cloud_models>` are set. 

After that, the radiative transfer scheme is set. There are currently two 
different options:

  - ``scm`` - the short characteristic method, available for CPU and GPU

  - ``disort`` - the discrete-ordinate solver DISORT, only available for CPU

Finally, the different :ref:`chemistry models <sec:chemistry_models>`, chemical species, and the
:ref:`opacity sources <sec:opacity_data>` are selected. It is important to note that
chemical species that are used as part of the chemistry models, should normally also have
an associated opacity source. Otherwise, the impact of that species on the 
resulting spectrum might be negligible and, therefore, its abundance will be 
likely be unconstrained.

On the other hand, it is theoretically also possible to add opacity species
without including this species in any of the chemistry models. In this case, the
abundance of this species will be zero and, thus, won't show up in the spectrum.


Model postprocessing
....................

The optional postprocessing file ``post_process.config`` has the following structure and
options:

.. include:: ../examples/post_process_emission.config
   :literal:

The default options that are used when this file is not present are:

  - ``Delete unused MultiNest files`` : No

  - ``Compute spectra`` : Yes

  - ``Save temperature structures`` : Yes

  - ``Save effective temperatures`` : Yes

  - ``Save contribution functions`` : No

  - ``Save chemical species`` : None

