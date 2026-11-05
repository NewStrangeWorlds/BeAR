
.. _sec:cloud_models:

Cloud Models
============

BeAR currently includes the following cloud descriptions:

  - :ref:`Grey cloud layer <sec:cloud_model_grey>`

  - :ref:`Non-grey cloud fit <sec:cloud_model_kh>` from 
    `Kitzmann \& Heng (2018) <https://ui.adsabs.harvard.edu/abs/2018MNRAS.475...94K/>`_

  - :ref:`Power law cloud layer <sec:cloud_model_power_law>` 

The code also supports the :ref:`mixing of different cloud models <sec:cloud_model_mixing>` in a single retrieval.

Cloud models are configured through the ``clouds`` array in the ``forward_model.toml`` file. To run a
retrieval without any clouds, set an empty array: ``clouds = []``.


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
