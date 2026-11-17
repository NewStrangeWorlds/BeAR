
.. _sec:cloud_models:

Cloud Models
============

BeAR currently includes the following cloud descriptions:

  - :ref:`Grey cloud layer <sec:cloud_model_grey>`

  - :ref:`Non-grey cloud fit <sec:cloud_model_kh>` from 
    `Kitzmann \& Heng (2018) <https://ui.adsabs.harvard.edu/abs/2018MNRAS.475...94K/>`_

  - :ref:`Power law cloud layer <sec:cloud_model_power_law>`

  - :ref:`Tabulated Mie cloud model <sec:cloud_model_mie>` with scattering properties from a
    pre-computed table of Mie calculations

The code also supports the :ref:`mixing of different cloud models <sec:cloud_model_mixing>` in a single retrieval.

Cloud models are configured through the ``clouds`` array in the ``forward_model.toml`` file. To run a
retrieval without any clouds, set an empty array: ``clouds = []``.

Each cloud model provides the optical depth, the single scattering albedo, and the asymmetry parameter of
the cloud for every atmospheric layer and spectral point. The grey, non-grey, and power-law models describe
purely absorbing clouds, unless the :ref:`scattering option <sec:cloud_model_grey_scattering>` of the grey
model is used, while the ``mie`` model always provides a scattering cloud. The scattering properties are
only used by the emission, secondary eclipse, and phase-curve forward models in combination with a
:ref:`radiative transfer scheme <sec:radiative_transfer>` that includes scattering. The default short
characteristics scheme treats clouds as pure absorbers. The transmission spectrum model uses the total cloud
optical depth as extinction, which is the appropriate treatment for the transit geometry.


.. _sec:cloud_model_grey:

Grey Cloud Model
----------------

The grey cloud layer model is the most simple cloud parametrisation available in BeAR.
It puts a cloud layer with a certain geometric height and a wavelength-independent
optical depth into the atmosphere. Cloud models are configured through the ``clouds`` array
in the ``forward_model.toml`` file. The grey model is selected by using the model keyword
:code:`grey`:

.. code::

   clouds = [ { model = "grey", params = [] } ]

The location of the cloud layer is set by a cloud top pressure :math:`p_\mathrm{top}`
that is one of the free parameters of the model. The cloud bottom pressure
:math:`p_\mathrm{bot}` is determined by a free parameter :math:`b`, such that
:math:`p_\mathrm{bot} = p_\mathrm{top} \, b`. The optical depth is divided between all
atmospheric layers that lie between these two pressures.

A grey cloud model has in total three free parameters that have to be added to the prior config
file:

  - (vertical) grey optical depth (``cloud_optical_depth``)

  - cloud top pressure in bar (``cloud_top_pressure``)

  - cloud bottom parameter :math:`b` (``cloud_bottom_fraction``)

The cloud priors in ``priors.config`` are name-keyed (via the names listed in parentheses above)
and therefore order-independent.

It is sometimes not possible to constrain the bottom of the cloud deck properly.
This is often the case when performing retrieval calculations for a transmission spectrum.
In this case, BeAR has the option to use a fixed cloud bottom. This can be enabled by
adding the :code:`fb` string to the ``params`` of the cloud model:

.. code::

   clouds = [ { model = "grey", params = ["fb"] } ]

By choosing this option, BeAR fixes the cloud bottom pressure at one atmospheric scale height
below the cloud top. In this case, only two free parameters are required in the prior file:

  - (vertical) grey optical depth (``cloud_optical_depth``)

  - cloud top pressure in bar (``cloud_top_pressure``)


.. _sec:cloud_model_grey_scattering:

Scattering grey cloud
.....................

By default, the grey cloud is a pure absorber. With the :code:`scattering` option, the cloud layer
scatters with a constant, wavelength-independent single scattering albedo :math:`\omega` and a
Henyey-Greenstein asymmetry parameter :math:`g`, both of which become additional free parameters:

.. code::

   clouds = [ { model = "grey", params = ["scattering"] } ]

The option can be combined with the fixed cloud bottom, e.g. ``params = ["fb", "scattering"]``. With the
``scattering`` option, the following two free parameters are required in addition to the ones listed above:

  - single scattering albedo :math:`\omega` (``cloud_single_scattering_albedo``), between 0 and 1

  - asymmetry parameter :math:`g` (``cloud_asym_param``), between -1 and 1

Values outside of these ranges are clipped internally. The scattering properties only have an effect in
combination with a radiative transfer scheme that includes scattering, see :ref:`sec:radiative_transfer`.
With the short characteristics scheme or in the transmission model, the cloud still acts as pure extinction.


.. _sec:cloud_model_kh:

Kitzmann & Heng Non-Grey Cloud Model
------------------------------------

This description of a non-grey cloud layer is based on analytical fits to the Mie efficiencies
of the cloud particles presented by 
`Kitzmann \& Heng (2018) <https://ui.adsabs.harvard.edu/abs/2018MNRAS.475...94K/>`_.
The Mie extinction efficiencies as a function of wavelength :math:`Q_\mathrm{ext}` 
are described via:

.. math:: Q_\mathrm{ext}(\lambda) = \frac{Q_1}{Q_0 x^{-a_0}_\lambda + x^{0.2}_\lambda} ,

where :math:`x_\lambda` is the wavelength-dependent size parameter of the cloud particles, 
:math:`Q_1` is a normalisation constant, :math:`Q_0` determines the :math:`x_\lambda` value 
at which :math:`Q_\mathrm{ext}` is peaking, and :math:`a_0` is the powerlaw index in the 
small particle limit, where Mie theory converges to the limit of Rayleigh scattering. 
The size parameter is a function of the particle radius :math:`a` and the wavelength:

.. math:: x_\lambda = 2 \pi a / \lambda .

The optical depth :math:`\tau_\lambda` of the cloud is determined based on the 
number density of cloud particles :math:`n_c` and the vertical extent of the cloud
:math:`\Delta z`:

.. math:: \tau(\lambda) = Q_\mathrm{ext}(\lambda) \pi a^2 n_c \Delta z .

Since, however, it is usually difficult to put good priors on the number density 
:math:`n_c`, BeAR instead uses the vertical optical depth at a reference wavelength
:math:`\lambda_\mathrm{ref}` to describe the optical depths as a function of wavelength:

.. math:: \tau(\lambda) = \tau(\lambda_\mathrm{ref}) 
   \frac{Q_\mathrm{ext}(\lambda)}{Q_\mathrm{ext}(\lambda_\mathrm{ref})}
   = \tau(\lambda_\mathrm{ref}) 
     \frac{Q_0 x^{-a_0}_{\lambda_{\mathrm{ref}}} + x^{0.2}_{\lambda_\mathrm{ref}}}{Q_0 x^{-a_0}_\lambda + x^{0.2}_\lambda}  .

By using the reference optical depth, the normalisation constant :math:`Q_1` is not required anymore.

The non-grey cloud description is selected by using the model keyword
:code:`KHnongrey`, with the reference wavelength in :math:`\mu\mathrm{m}` given as the first
entry of ``params`` in the ``clouds`` array in the ``forward_model.toml`` file:

.. code::

   clouds = [ { model = "KHnongrey", params = ["1.0"] } ]

The reference wavelength does not need to be within the wavelength range of the
retrieval calculations.

A non-grey cloud model has in total six free parameters that have to be added to the prior config
file:

  - (vertical) optical depth at the reference wavelength (``cloud_optical_depth``)

  - :math:`Q_0` (``cloud_q0``)

  - :math:`a_0` (``cloud_a0``)

  - particle size :math:`a` in microns (``cloud_particle_size``)

  - cloud top pressure in bar (``cloud_top_pressure``)

  - cloud bottom parameter :math:`b` (``cloud_bottom_fraction``)

The cloud priors in ``priors.config`` are name-keyed (via the names listed in parentheses above)
and therefore order-independent.

The location of the cloud layer is set by a cloud top pressure :math:`p_\mathrm{top}`
that is one of the free parameters of the model. The cloud bottom pressure
:math:`p_\mathrm{bot}` is determined by a free parameter :math:`b`, such that
:math:`p_\mathrm{bot} = p_\mathrm{top} \, b`. The optical depth is divided between all
atmospheric layers that lie between these two pressures.

It is sometimes not possible to constrain the bottom of the cloud deck properly.
This is often the case when performing retrieval calculations for a transmission spectrum.
In this case, BeAR has the option to use a fixed cloud bottom. This can be enabled by
adding the :code:`fb` string to ``params`` after the reference wavelength:

.. code::

   clouds = [ { model = "KHnongrey", params = ["1.0", "fb"] } ]

By choosing this option, BeAR fixes the cloud bottom pressure at one atmospheric scale height
below the cloud top. In this case, only five free parameters are required in the prior file:

  - (vertical) optical depth at the reference wavelength (``cloud_optical_depth``)

  - :math:`Q_0` (``cloud_q0``)

  - :math:`a_0` (``cloud_a0``)

  - particle size :math:`a` (``cloud_particle_size``)

  - cloud top pressure in bar (``cloud_top_pressure``)


.. _sec:cloud_model_power_law:

Power-Law Cloud Model
---------------------

This cloud model uses a power law to describe the wavelength-dependent optical depth of 
the cloud layer. It uses the optical depth at a reference wavelength as normalisation:

.. math:: \tau(\lambda) = \tau(\lambda_\mathrm{ref}) \frac{\lambda^e}{\lambda_\mathrm{ref}^e} ,

where :math:`e` is the exponent of the power law. To simulate Rayleigh scattering, for
example, :math:`e=-4`.

The power-law cloud model is selected by using the model keyword
:code:`power_law`, with the reference wavelength in :math:`\mu\mathrm{m}` given as the first
entry of ``params`` in the ``clouds`` array in the ``forward_model.toml`` file:

.. code::

   clouds = [ { model = "power_law", params = ["1.0"] } ]

The reference wavelength does not need to be within the wavelength range of the
retrieval calculations.

A power-law cloud model has in total four free parameters that have to be added to the prior config
file:

  - (vertical) optical depth at the reference wavelength (``cloud_optical_depth``)

  - power-law exponent :math:`e` (``cloud_exponent``)

  - cloud top pressure in bar (``cloud_top_pressure``)

  - cloud bottom parameter :math:`b` (``cloud_bottom_fraction``)

The cloud priors in ``priors.config`` are name-keyed (via the names listed in parentheses above)
and therefore order-independent.

The location of the cloud layer is set by a cloud top pressure :math:`p_\mathrm{top}`
that is one of the free parameters of the model. The cloud bottom pressure
:math:`p_\mathrm{bot}` is determined by a free parameter :math:`b`, such that
:math:`p_\mathrm{bot} = p_\mathrm{top} \, b`. The optical depth is divided between all
atmospheric layers that lie between these two pressures.

It is sometimes not possible to constrain the bottom of the cloud deck properly.
This is often the case when performing retrieval calculations for a transmission spectrum.
In this case, BeAR has the option to use a fixed cloud bottom. This can be enabled by
adding the :code:`fb` string to ``params`` after the reference wavelength:

.. code::

   clouds = [ { model = "power_law", params = ["1.0", "fb"] } ]

By choosing this option, BeAR fixes the cloud bottom pressure at one atmospheric scale height
below the cloud top. In this case, only three parameters are required in the prior file:

  - (vertical) optical depth at the reference wavelength (``cloud_optical_depth``)

  - power-law exponent :math:`e` (``cloud_exponent``)

  - cloud top pressure in bar (``cloud_top_pressure``)


.. _sec:cloud_model_mie:

Tabulated Mie Cloud Model
-------------------------

This cloud model describes a non-grey cloud layer whose optical properties are taken from a
pre-computed table of Mie calculations for a given cloud particle material and log-normal particle
size distributions. The table provides the mean extinction cross section per particle
:math:`C_\mathrm{ext}`, the single scattering albedo :math:`\omega`, and the asymmetry parameter
:math:`g` as a function of the median particle radius :math:`r_\mathrm{g}`, the geometric standard
deviation :math:`\sigma_\mathrm{g}` of the size distribution, and the wavelength. In contrast to the
other cloud models, this model therefore yields a scattering cloud with a wavelength-dependent single
scattering albedo and asymmetry parameter.

Like the non-grey and power-law models, the optical depth is parametrised by its value at a reference
wavelength :math:`\lambda_\mathrm{ref}`:

.. math:: \tau(\lambda) = \tau(\lambda_\mathrm{ref})
   \frac{C_\mathrm{ext}(\lambda)}{C_\mathrm{ext}(\lambda_\mathrm{ref})} ,

so that the number density of the cloud particles is not required as a parameter. The single scattering
albedo and the asymmetry parameter are taken directly from the table. The table is interpolated
bilinearly in :math:`\ln r_\mathrm{g}` and :math:`\sigma_\mathrm{g}` and linearly in
:math:`\ln \lambda`, with the extinction cross section interpolated in logarithmic space. Particle sizes
and distribution widths outside of the tabulated range are clamped to the table boundaries. Spectral
grid points outside of the wavelength range of the table use the values at the table edges; BeAR prints
a warning in this case.

The model is selected with the keyword :code:`mie`. The first entry of ``params`` is the file name of
the table, relative to the retrieval folder, and the second one the reference wavelength in
:math:`\mu\mathrm{m}`, which has to lie inside the wavelength range of the table:

.. code::

   clouds = [ { model = "mie", params = ["mg2sio4_mie.dat", "1.0"] } ]

Two optional keywords can be appended to ``params``:

  - :code:`fb` - fixes the cloud bottom one atmospheric scale height below the cloud top, as for the
    other cloud models

  - :code:`sigma_g=<value>` - fixes the width of the size distribution to the given value instead of
    retrieving it, for example ``"sigma_g=1.5"``. The value has to lie inside the range of the table.

For example:

.. code::

   clouds = [ { model = "mie", params = ["mg2sio4_mie.dat", "1.0", "fb", "sigma_g=1.5"] } ]

The Mie cloud model has up to five free parameters that have to be added to the prior config file:

  - (vertical) optical depth at the reference wavelength (``cloud_optical_depth``)

  - median particle radius :math:`r_\mathrm{g}` in microns (``cloud_particle_size``)

  - geometric standard deviation :math:`\sigma_\mathrm{g}` of the size distribution (``cloud_sigma_g``),
    only if it is not fixed with the ``sigma_g=`` keyword

  - cloud top pressure in bar (``cloud_top_pressure``)

  - cloud bottom parameter :math:`b` (``cloud_bottom_fraction``), only without the ``fb`` keyword

The location and the vertical extent of the cloud are treated in the same way as for the grey cloud
model: the optical depth at the reference wavelength is distributed evenly over all atmospheric layers
between the cloud top and bottom.

The scattering properties of this cloud model only have an effect in combination with a radiative
transfer scheme that includes scattering, see :ref:`sec:radiative_transfer`.


Table format
............

The table is a plain-text file in which lines starting with ``#`` are treated as comments and can
appear anywhere. It has the following layout:

.. code::

   nb_r nb_sigma nb_wavelength
   r_g values          (one line, in micron, ascending)
   sigma_g values      (one line, ascending)
   wavelength values   (one line, in micron, ascending)
   for each r_g (outer loop) and each sigma_g (inner loop):
     C_ext values      (one line, in cm^2, one value per wavelength)
     omega values      (one line, one value per wavelength)
     g values          (one line, one value per wavelength)

The three axes have to be strictly ascending and the extinction cross sections have to be positive. The
Python script ``tools/mie_table_to_ascii.py`` in the BeAR repository converts an HDF5 file with the
datasets ``r_g_micron``, ``sigma_g``, ``wavelength_micron``, ``c_ext_cm2``, ``single_scattering_albedo``,
and ``asymmetry_parameter`` into this format. The wavelength range can optionally be cropped to reduce
the file size:

.. code:: bash

   python3 tools/mie_table_to_ascii.py mg2sio4_optics_grid.h5 mg2sio4_mie.dat --wavelength-min 0.3 --wavelength-max 30


.. _sec:cloud_model_mixing:

Mixing different cloud models
-----------------------------

BeAR also has the ability to use multiple clouds simultaneously. For example, to perform a retrieval with
grey cloud layer and a second cloud with a power-law optical depth to simulate, for example, a haze-like
behaviour, several entries can be added to the ``clouds`` array in the ``forward_model.toml`` file:

.. code::

   clouds = [ { model = "grey", params = [] }, { model = "power_law", params = ["1.0", "fb"] } ]

BeAR will call the cloud models in the order they appear in this array. Clouds can also be overlapping
in pressure.

The cloud priors in ``priors.config`` are name-keyed and order-independent. When more than one
cloud model is used, the parameter names of each cloud are suffixed with the 1-based index of the
cloud in the ``clouds`` array (``_1``, ``_2``, ...) so that they stay unique. With a single cloud
model, the names are used without a suffix. For the example above, the following priors need to be
listed:

  - grey optical depth (``cloud_optical_depth_1``)

  - cloud top pressure in bar (``cloud_top_pressure_1``)

  - cloud bottom parameter (``cloud_bottom_fraction_1``)

  - optical depth at the reference wavelength (``cloud_optical_depth_2``)

  - power-law exponent (``cloud_exponent_2``)

  - cloud top pressure in bar (``cloud_top_pressure_2``)
