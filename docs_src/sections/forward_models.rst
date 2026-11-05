
.. _sec:forward_models:

Forward Models
==============

BeAR currently includes the following forward models:

  - :ref:`Transmission spectrum <sec:forward_model_transmission>`

  - :ref:`Emission spectrum <sec:forward_model_em>`

  - :ref:`Secondary eclipse spectrum <sec:forward_model_se>`

  - :ref:`Secondary eclipse spectra with a planetary blackbody <sec:forward_model_se_bb>`

  - :ref:`Flat line <sec:forward_model_flat>`

  - Phase curve - the high-resolution phase-curve model, documented in a later section.

The secondary eclipse blackbody and flat line models are usually used to check if the
observational data is more likely explained by simpler models, such as a flat line.

The forward model is selected in the main retrieval config file ``retrieval.toml`` via the
``forward_model_type`` key in the ``[retrieval]`` table:

.. code:: toml

   [retrieval]
   forward_model_type = "emission"

The accepted values (each with a short alias) are:

  - ``emission`` / ``em``

  - ``secondary_eclipse`` / ``se``

  - ``secondary_eclipse_bb`` / ``se_bb``

  - ``transmission`` / ``trans``

  - ``flat_line`` / ``fl``

  - ``phase_curve`` / ``pc``

Each model (except the flat line) requires its own config file ``forward_model.toml`` that
is discussed below.


Model postprocessing
--------------------

After the retrieval calculations are finished, BeAR will perform a postprocessing on
the resulting posterior sample. Depending on the chosen forward model, different
postprocessing steps can be used, including saving all posterior spectra or temperature
profiles. Details can be found in the :ref:`section <sec:output_files>`
on output files.

These postprocessing steps can also be configured in the optional
``post_process.toml`` file. If this file is not present, BeAR will use
default settings for the postprocessing. Depending on the forward model,
available options are currently:

  - ``delete_sampler_files`` - This option will delete all MultiNest sampler files that
    are not used in the post processing. Boolean value: ``true`` or ``false``.

  - ``save_spectra`` - Compute the spectra for each model in the posterior sample. Spectra
    will be saved for each observation/instrument individually. Additionally, a high-resolution
    spectrum of the best-fit model will be computed and saved. Boolean value: ``true`` or ``false``.

  - ``save_temperatures`` - Compute and save the temperature profile for each model in the posterior
    sample. Boolean value: ``true`` or ``false``.

  - ``save_effective_temperatures`` - Compute and save the effective temperature for each model in the
    posterior sample. Boolean value: ``true`` or ``false``.

  - ``save_contribution_functions`` - Compute and save the contribution functions for each
    observation/instrument for the best-fit model. Boolean value: ``true`` or ``false``.

  - ``species_to_save`` - Save the mixing ratio profiles for selected chemical species for all
    posterior samples. The value is an array of the formulas of the chemical species that should be
    saved. For example, ``["H2O", "CO2"]`` will save the mixing ratios of water and carbon dioxide.
    Species that BeAR does not know will be ignored.



.. _sec:forward_model_flat:

Flat line
---------

This model simply fits a flat line through the observational data. Usually, this
model is used to test if the use of a more complex model is warranted to explain
the observation. The test is normally done by comparing their
Bayesian evidences.

In the main retrieval config file ``retrieval.toml`` it is selected by choosing:

.. code:: toml

   [retrieval]
   forward_model_type = "flat_line"


Model config file
.................

The flat line model does not need a ``forward_model.toml`` file since it
has no configurable parameters. In the prior distribution file this model
requires a single free parameter, ``spectrum_value``, that determines the flat line.

This parameter needs to have the same units as the observational data. For example,
in case of a transmission spectrum, this parameter refers to the transit depth in ppm, while
for an emission spectrum, it needs to have units of
:math:`\mathrm{W} \mathrm{m^{-2}} \mathrm{\mu m^{-1}}`.


Model postprocessing
....................

The optional postprocessing file ``post_process.toml`` has the following structure and
options:

.. include:: ../examples/post_process_flatline.toml
   :literal:

The default options that are used when this file is not present are:

  - ``delete_sampler_files`` : false

  - ``save_spectra`` : false



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

In the retrieval config file ``retrieval.toml`` it is selected by choosing:

.. code:: toml

   [retrieval]
   forward_model_type = "transmission"

The transmission spectrum forward model has three general free parameters that
have to be added to the priors configuration file. Since ``priors.config`` is
name-keyed, the entries are identified by their names rather than their order:

  - ``log_g`` - logarithm of surface gravity :math:`\log g` in cgs units

  - ``bottom_radius`` - planet radius

  - ``star_radius`` - stellar radius


Model config file
.................

The ``forward_model.toml`` file for the transmission spectrum model has
the following structure:

.. include:: ../examples/forward_model_transmission.toml
   :literal:

The file starts with two inline keys. The first, ``fit_mode``, selects the retrieval of the
mean molecular weight or the scale height. Usually, BeAR will determine the mean molecular
weight based on the chemical abundances of the species included in the retrieval.
Sometimes, however, the background species in an atmosphere might not be
known. For such a case, BeAR can use the mean molecular weight as a free
parameter.

Furthermore, for transmission spectra, the surface gravity, the mean molecular
weight, and the temperature might become degenerate in a retrieval. This is
often the case when no constraints on the surface gravity can be provided
and the dominating background species is not known. In such a scenario, BeAR
can use the (constant) atmospheric scale height in units of km as a free
parameter.

For the standard case, ``fit_mode = "no"`` is used. If BeAR should
use the mean molecular weight as a free parameter, ``fit_mode = "mmw"`` is used, while
``fit_mode = "sh"`` is used when the scale height should be used as a free parameter instead.
If either the mean molecular weight or the scale height is chosen as a free parameter, a
corresponding prior (``mean_molecular_weight`` or ``scale_height``) needs to be added to the
prior configuration file.

The second key, ``use_variable_gravity``, is a boolean (``true`` or ``false``) that determines
whether the gravitational acceleration is kept constant or allowed to vary with altitude
throughout the atmosphere.

The ``atmosphere`` table sets the vertical discretisation of the atmosphere. This includes the
number of atmospheric layers (``nb_grid_points``) as well as the bottom (``bottom_pressure``) and
top-of-atmosphere (``top_pressure``) pressures in units of bars.

With the next two entries, the parametrisation for the
:ref:`temperature profile <sec:temperature_profiles>` and
the :ref:`cloud models <sec:cloud_models>` are set. After that, optional
modules can be added to the forward model via the ``modules`` array.

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

The optional postprocessing file ``post_process.toml`` has the following structure and
options:

.. include:: ../examples/post_process_transmission.toml
   :literal:

The default options that are used when this file is not present are:

  - ``delete_sampler_files`` : false

  - ``save_spectra`` : true

  - ``save_temperatures`` : false

  - ``species_to_save`` : none



.. _sec:forward_model_se:

Secondary eclipse spectrum
---------------------------

This forward model computes the wavelength-dependent secondary eclipse (or occultation)
depth :math:`D(\lambda)` of an exoplanet atmosphere, given by

.. math::
  D(\lambda) = \frac{F_p(\lambda)}{F_*(\lambda)} \left(\frac{R_p}{R_*}\right)^2 \ ,

where :math:`F_p(\lambda)` is the outgoing flux at top of the planet's atmosphere,
:math:`F_*(\lambda)` is the stellar photospheric flux, :math:`R_p`
is the planetary radius and :math:`R_*` the radius of the host star.
In BeAR, :math:`D(\lambda)` has units of ppm.

In the retrieval config file ``retrieval.toml`` it is selected by choosing:

.. code:: toml

   [retrieval]
   forward_model_type = "secondary_eclipse"

The secondary eclipse spectrum forward model has two general free parameters that
have to be added to the priors configuration file. Since ``priors.config`` is
name-keyed, the entries are identified by their names:

  - ``log_g`` - logarithm of surface gravity :math:`\log g` in cgs units

  - ``radius_ratio`` - ratio of the planet's and stellar radius :math:`\mathrm{R_p/R_*}`


Model config file
.................

The ``forward_model.toml`` file for the secondary eclipse spectrum model has
the following structure:

.. include:: ../examples/forward_model_se.toml
   :literal:

The ``atmosphere`` table sets the vertical discretisation of the atmosphere. This includes the
number of atmospheric layers (``nb_grid_points``) as well as the bottom (``bottom_pressure``) and
top-of-atmosphere (``top_pressure``) pressures in units of bars.

With the next entries, the parametrisation for the :ref:`temperature profile <sec:temperature_profiles>`,
the :ref:`stellar spectrum <sec:stellar_spectra>`, and the :ref:`cloud models <sec:cloud_models>` are set.

After that, the radiative transfer scheme is chosen via
``radiative_transfer = { model = "...", params = [...] }``. There are currently three
different options:

  - ``scm`` - the short characteristic method (ShortCharacteristics), ``params = []``,
    available for CPU and GPU.

  - ``disort`` - the discrete-ordinate solver DisORT++ (an external library), which takes
    exactly one parameter, the number of streams (must be 4, 8, 16, or 32). Only available
    for CPU; it throws an error if run on the GPU.

  - ``adding_doubling`` - the adding-doubling method (AddingDoublingRT), which takes exactly
    one parameter, the number of quadrature points. Available for CPU and GPU, and supports
    aerosol/multiple scattering.

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

The optional postprocessing file ``post_process.toml`` has the following structure and
options:

.. include:: ../examples/post_process_se.toml
   :literal:

The default options that are used when this file is not present are:

  - ``delete_sampler_files`` : false

  - ``save_spectra`` : true

  - ``save_temperatures`` : true

  - ``save_contribution_functions`` : false

  - ``species_to_save`` : none



.. _sec:forward_model_se_bb:

Secondary eclipse spectrum with planetary blackbody
---------------------------------------------------

This forward model computes the wavelength-dependent secondary eclipse (or occultation)
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

In the retrieval config file ``retrieval.toml`` it is selected by choosing:

.. code:: toml

   [retrieval]
   forward_model_type = "secondary_eclipse_bb"

This forward model has two general free parameters that
have to be added to the priors configuration file:

  - ``planet_temperature`` - the planet's effective temperature in Kelvin

  - ``radius_ratio`` - ratio of the planet's and stellar radius :math:`\mathrm{R_p/R_*}`

In addition, the parameters of the chosen stellar spectrum model are required, e.g.
``star_t_eff`` for a blackbody star (see :ref:`sec:stellar_spectra`).


Model config file
.................

The ``forward_model.toml`` file for the secondary eclipse spectrum model has
the following structure:

.. include:: ../examples/forward_model_se_bb.toml
   :literal:

It contains a single option related to the description of the stellar spectrum. Information
on the available options can be found in the :ref:`section <sec:stellar_spectra>` on stellar spectra.


Model postprocessing
....................

The optional postprocessing file ``post_process.toml`` has the following structure and
options:

.. include:: ../examples/post_process_sebb.toml
   :literal:

The default options that are used when this file is not present are:

  - ``delete_sampler_files`` : false

  - ``save_spectra`` : true



.. _sec:forward_model_em:

Emission spectrum
-----------------

This forward model computes the emission spectrum :math:`F(\lambda)` of an exoplanet or brown dwarf atmosphere.
In BeAR, :math:`F(\lambda)` has units of :math:`\mathrm{W} \mathrm{m^{-2}} \mathrm{\mu m^{-1}}`.

In the retrieval config file ``retrieval.toml`` it is selected by choosing:

.. code:: toml

   [retrieval]
   forward_model_type = "emission"

The emission spectrum forward model has three general free parameters that
have to be added to the priors configuration file. Since ``priors.config`` is
name-keyed, the entries are identified by their names:

  - ``log_g`` - logarithm of surface gravity :math:`\log g` in cgs units

  - ``scaling_factor`` - a scaling factor :math:`f` for the radius/distance relationship,
    where the radius is internally set to 1 Jupiter radius

  - ``distance`` - distance to the object


Model config file
.................

The ``forward_model.toml`` file for the emission spectrum model has
the following structure:

.. include:: ../examples/forward_model_em.toml
   :literal:

The ``atmosphere`` table sets the vertical discretisation of the atmosphere. This includes the
number of atmospheric layers (``nb_grid_points``) as well as the bottom (``bottom_pressure``) and
top-of-atmosphere (``top_pressure``) pressures in units of bars.

With the next two entries, the parametrisation for the
:ref:`temperature profile <sec:temperature_profiles>` and
the :ref:`cloud models <sec:cloud_models>` are set.

After that, the radiative transfer scheme is set via
``radiative_transfer = { model = "...", params = [...] }``. There are currently three
different options:

  - ``scm`` - the short characteristic method (ShortCharacteristics), ``params = []``,
    available for CPU and GPU.

  - ``disort`` - the discrete-ordinate solver DisORT++ (an external library), which takes
    exactly one parameter, the number of streams (must be 4, 8, 16, or 32). Only available
    for CPU; it throws an error if run on the GPU.

  - ``adding_doubling`` - the adding-doubling method (AddingDoublingRT), which takes exactly
    one parameter, the number of quadrature points. Available for CPU and GPU, and supports
    aerosol/multiple scattering.

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

The optional postprocessing file ``post_process.toml`` has the following structure and
options:

.. include:: ../examples/post_process_emission.toml
   :literal:

The default options that are used when this file is not present are:

  - ``delete_sampler_files`` : false

  - ``save_spectra`` : true

  - ``save_temperatures`` : true

  - ``save_effective_temperatures`` : true

  - ``save_contribution_functions`` : false

  - ``species_to_save`` : none


.. _sec:forward_model_phase_curve:

Phase curve (high-resolution)
-----------------------------

This forward model is used for high-resolution cross-correlation spectroscopy and
phase-curve retrievals. It computes both a low-resolution occultation depth and a
high-resolution planet-to-star flux ratio :math:`F_\mathrm{p}/F_\mathrm{s}` on a
separate high-resolution spectral grid. The full high-resolution workflow — the
orbital-velocity model, exposure blurring, model filtering and re-injection, and
the cross-correlation likelihoods — is described on a dedicated page:
:ref:`sec:high_resolution`.

In the retrieval config file ``retrieval.toml`` it is selected by choosing:

.. code:: toml

   [retrieval]
   forward_model_type = "phase_curve"

It has two general free parameters, identified by name in ``priors.config``:

  - ``log_g`` - logarithm of surface gravity :math:`\log g` in cgs units

  - ``radius_ratio`` - the planet-to-star radius ratio :math:`R_\mathrm{p}/R_\mathrm{s}`

When high-resolution observations are present, BeAR additionally appends the
high-resolution tail parameters ``kp``, ``vsys``, ``dphi`` and, optionally,
``alpha`` to the parameter list, together with any parameters required by the
selected modules (see :ref:`sec:high_resolution` for details).


Model config file
.................

The ``forward_model.toml`` file for the phase-curve model uses the same building
blocks as the other models (``atmosphere``, ``temperature``, ``stellar_spectrum``,
``clouds``, ``radiative_transfer``, ``chemistry``, ``opacity``) and adds two
high-resolution-specific pieces:

  - an optional ``opacity_highres`` block, a separate opacity table for the
    high-resolution grid (see :ref:`sec:opacity_highres`); if omitted it falls
    back to the normal ``opacity`` list, and

  - an optional ``modules`` array for post-processing effects such as velocity
    broadening, phase-resolved broadening, or stellar contamination
    (see :ref:`sec:high_resolution`).

The ``stellar_spectrum`` table of this model additionally accepts the
``highres_stellar_smooth_sigma`` key to smooth the stellar spectrum for the
high-resolution model (see :ref:`sec:stellar_spectra`). Its ``post_process.toml``
file supports the same keys as the secondary eclipse model plus the boolean
``save_hr_spectra_all`` (default ``false``), which writes the high-resolution
spectrum for every posterior sample instead of only the best-fit model.

The ``stellar_spectrum`` table may also carry the optional
``highres_stellar_smooth_sigma`` key used for high-resolution re-injection (see
:ref:`sec:high_resolution`).
