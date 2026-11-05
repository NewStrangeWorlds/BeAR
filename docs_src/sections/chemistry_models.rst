
.. _sec:chemistry_models:

Chemistry Models
=================

BeAR currently includes the following chemistry models:

  - :ref:`Equilibrium chemistry <sec:chemistry_model_eq>` (:code:`eq` / :code:`equilibrium`)

  - :ref:`Background species <sec:chemistry_model_bg>` (:code:`bg` / :code:`background`)

  - :ref:`Isoprofiles <sec:chemistry_model_iso>` (:code:`iso` / :code:`isoprofile`)

  - :ref:`Isoprofiles with centred-log-ratio priors <sec:chemistry_model_iso_clr>` (:code:`iso_clr` / :code:`isoprofile_clr`)

  - :ref:`Isoprofiles with a retrieved isotopologue ratio <sec:chemistry_model_iso_ratio>` (:code:`iso_ratio` / :code:`isoprofile_ratio`)

  - :ref:`Piecewise polynomials <sec:chemistry_model_poly>` (:code:`free`)

  - :ref:`Cubic b splines <sec:chemistry_model_cs>` (:code:`free_cs` / :code:`free_cspline`)

  - :ref:`Piecewise cubic Hermite interpolation <sec:chemistry_model_pchip>` (:code:`free_p` / :code:`free_pchip`)

  - :ref:`Step function <sec:chemistry_model_sf>` (:code:`sf` / :code:`step_function`)

Details on the implementation of these parametrisations are
discussed :ref:`here <sec:profile_parametrisations>`.

BeAR can also mix several different of these models in a single retrieval. This
is explained :ref:`in this section <sec:chemistry_model_mixing>`. All chemical
species that BeAR uses have to be defined in a corresponding header file.
This is described at the :ref:`bottom <sec:chemical_species>` of this chapter.


.. _sec:chemistry_model_iso:

Isoprofiles
-----------

In the corresponding ``forward_model.toml`` of the chosen forward model,
the isoprofile chemistry model is chosen by using the keyword :code:`iso`
(or :code:`isoprofile`), followed by a list of species that should be used by the retrieval:

.. code:: toml

   chemistry = [
     { model = "iso", params = ["H2O", "CH4", "NH3", "K", "H2S", "CO2"] },
   ]

For each species, a corresponding entry in the prior config file is required.
The priors are name-keyed and therefore order-independent: for the isoprofile model,
each species uses the prior name :code:`chem_<species>_mr`, with the species symbol written
in lower case. For the example above, the priors :code:`chem_h2o_mr`, :code:`chem_ch4_mr`,
:code:`chem_nh3_mr`, :code:`chem_k_mr`, :code:`chem_h2s_mr`, and :code:`chem_co2_mr` are expected.

Furthermore, as described by `Kitzmann et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020ApJ...890..174K/>`_,
the abundance of sodium (Na) is linked to the one of potassium (K) through the ratio of their
solar elemental abundances (:code:`solar_na_k = 16.2181`). Thus, in the example above, the Na abundance
would automatically be derived from K, without explicitly including the former as a free parameter.

If Na should be retrieved independently, it needs to be added as a separate species as shown below.

.. code:: toml

   chemistry = [
     { model = "iso", params = ["H2O", "CH4", "NH3", "K", "H2S", "CO2", "Na"] },
   ]


.. _sec:chemistry_model_iso_clr:

Isoprofiles with centred-log-ratio priors
-----------------------------------------

The standard isoprofile chemistry model described above fills up the background of the atmosphere with a mixture of molecular hydrogen and helium.
While this valid for gas giants or brown dwarfs, this is less suitable for atmospheres that are not dominated by H2 and He.
For cases, where the background gas is not explicitly known, centred-log-ratio priors for the chemical abundances can be used. These can more
efficiently sample the prior distributions in a way to determine which of the chosen species is the dominant one.
More detailed information on these priors can be found in `Benneke & Seager (2012) <https://ui.adsabs.harvard.edu/abs/2012ApJ...753..100B/>`_.

For a given mixture of :math:`n` gases, the centred-log-ratio conversion (clr) for the mixing ration :math:`x_j` of a given molecule :math:`j`
in the mixture is given by

.. math:: \xi_j = \mathrm{clr}(x_j) = \ln \frac{x_j}{g(\mathbf x)}

where :math:`g(\mathbf x)` is the geometric mean of all mixing ratios :math:`\mathbf{x}`:

.. math::
  g(\mathbf{x}) = \left( \prod_{j=1}^{n} x_j \right)^{1/n} \ .

Due to the constraint that

.. math::
  \sum_{j=1}^n x_j = 1 \quad \text{or}  \quad \sum_{j=1}^n \xi_j = 0 \ ,

only :math:`n-1` free parameters are needed in the retrieval.

The prior distributions for the free parameters should be uniform and vary between 0 and 1. BeAR will use these priors to
produce :math:`\xi_j` values subject to the constraints that :math:`\min\left(\mathbf{x}\right) = 10^{-12}` and :math:`\max\left(\mathbf{x}\right) = 1`.
It is important to note that the prior boundaries for :math:`\xi_j` depend on the number of molecules in the retrieval and the
chosen value of the smallest allowed mixing ratio. For uniform priors between 0 and 1, BeAR will do the conversion to :math:`\xi_j` and then
to the mixing ratios :math:`x_j` automatically.

This chemistry model is chosen in the ``forward_model.toml`` file by using

.. code:: toml

   chemistry = [
     { model = "iso_clr", params = ["H2O", "CH4", "NH3", "H2S", "CO2"] },
   ]

As noted above, for these 5 species only 4 priors are needed due to the constraint that the sum of their mixing ratios has to be unity.
This chemistry model will, therefore, select one of these free species as the dominant background gas. In the prior config file,
the clr parameters use the name :code:`chem_<species>_clr` (species symbol in lower case), and only priors for the first four molecules
would need to be listed. The abundance of the last species in the list above will be determined from the other molecules.


.. _sec:chemistry_model_iso_ratio:

Isoprofiles with a retrieved isotopologue ratio
-----------------------------------------------

This model is a variant of the isoprofile model in which the abundance of one
isotopologue is not retrieved directly but derived from its parent molecule through
a retrieved isotopic ratio. It is intended for high-resolution retrievals of, for
example, the :math:`^{13}\mathrm{CO}/^{12}\mathrm{CO}` ratio. The model is chosen
in the ``forward_model.toml`` file with the keyword :code:`iso_ratio`
(or :code:`isoprofile_ratio`), followed by the list of species:

.. code:: toml

   chemistry = [
     { model = "iso_ratio", params = ["H2O", "CO", "13C16O"] },
   ]

At least two species are required. The **last** species in the list is the isotopologue
and the **second-to-last** one is its parent molecule. All species except the last one
are treated as ordinary isoprofiles with constant mixing ratios, using the prior names
:code:`chem_<species>_mr` as in the isoprofile model. The last species uses the prior
name :code:`chem_<species>_ratio`. For the example above, the priors
:code:`chem_h2o_mr`, :code:`chem_co_mr`, and :code:`chem_13c16o_ratio` are expected.

The ratio parameter is the logarithmic offset from the terrestrial standard ratio,
following the convention of, e.g.,
`Line et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021Natur.598..580L/>`_:

.. math:: \left[ \mathrm{iso}/\mathrm{parent} \right] = \log_{10} \left( \frac{x_\mathrm{iso}}{x_\mathrm{parent}} \right) - \log_{10} \left( \frac{x_\mathrm{iso}}{x_\mathrm{parent}} \right)_\oplus \ ,

such that the mixing ratio of the isotopologue is given by

.. math:: x_\mathrm{iso} = \frac{x_\mathrm{parent}}{89} \, 10^{\left[ \mathrm{iso}/\mathrm{parent} \right]} \ .

A value of 0 therefore corresponds to the terrestrial ratio, while positive values
indicate an enrichment of the isotopologue. Note that the standard ratio of 89
(the terrestrial :math:`^{12}\mathrm{C}/^{13}\mathrm{C}` ratio) is currently
hard-coded in ``src/chemistry/iso_ratio_chemistry.h`` and is used regardless of
the chosen species.

In contrast to the standard isoprofile model, the Na abundance is *not* automatically
linked to K in this model. The remaining background of the atmosphere is treated in
the same way as for the other free chemistry models, see the
:ref:`section on mixing chemistry models <sec:chemistry_model_mixing>`.


.. _sec:chemistry_model_poly:

Piecewise polynomials
---------------------

If a chemical species, say, H2O should be included using non-constant abundances by using
the parametrisation with piecewise polynomials, the keyword :code:`free`, followed by the species'
formula, the number of elements :code:`k`, and the polynomial degree :code:`q` need to be added to
the ``forward_model.toml`` file. This model takes exactly three parameters:

.. code:: toml

   chemistry = [
     { model = "free", params = ["H2O", "k", "q"] },
   ]

Note, that unlike the isoprofile chemistry above, only a single species can be used here.
As stated :ref:`here <sec:profile_parametrisations>`, the number of free parameters for this
parametrisation is :math:`k q + 1`.
The free parameters correspond to the mixing ratios of the chosen species at the :math:`k q + 1`
discrete pressure points that are distributed equidistantly in logarithmic pressure space. In the
prior config file, they are name-keyed as :code:`chem_<species>_mr_<i>`, running from
:code:`chem_<species>_mr_0` at the bottom of the atmosphere to :code:`chem_<species>_mr_<kq>` at its top.


.. _sec:chemistry_model_cs:

Cubic b splines
---------------

For parametrising the vertical profile of a species, for example, CH4 by a cubic b spline, the
keyword :code:`free_cs` (or :code:`free_cspline`), followed by the species'
formula and the number of points :code:`k` need to be added to the ``forward_model.toml`` file. This
model takes exactly two parameters:

.. code:: toml

   chemistry = [
     { model = "free_cs", params = ["CH4", "k"] },
   ]

Note, that like for the piecewise polynomials, only a single species can be used here.
As stated :ref:`here <sec:profile_parametrisations>`, the number of free parameters for this
parametrisation is equal to the number of points :math:`k` and needs to be at least 5.
They correspond to the mixing ratios at the :math:`k` discrete pressure points that are distributed
equidistantly in logarithmic pressure space.

In the prior config file, these free parameters are name-keyed as :code:`chem_<species>_mr_<i>`,
running from :code:`chem_<species>_mr_0` at the bottom of the atmosphere to
:code:`chem_<species>_mr_<k-1>` at its top.


.. _sec:chemistry_model_pchip:

Piecewise cubic Hermite interpolation
-------------------------------------

For parametrising the vertical profile of a species, for example, CH4 by a monotone
piecewise cubic Hermite interpolant (PCHIP), the keyword :code:`free_p` (or :code:`free_pchip`),
followed by the species' formula and the number of control points :code:`k` need to be added
to the ``forward_model.toml`` file. This model takes exactly two parameters:

.. code:: toml

   chemistry = [
     { model = "free_p", params = ["CH4", "k"] },
   ]

Note, that like for the piecewise polynomials and the cubic b splines, only a single species
can be used here. The number of free parameters is equal to the number of control points
:math:`k`, which needs to be at least 4. They correspond to the mixing ratios at the :math:`k`
control points that are distributed equidistantly in logarithmic pressure space.

BeAR uses the monotone PCHIP interpolant from the Boost.Math library to construct the profile
between the control points. Unlike the :ref:`cubic b spline <sec:chemistry_model_cs>` model
(which requires at least 5 control points and does not preserve monotonicity), the PCHIP
interpolant is shape-preserving: it does not overshoot between adjacent control points and
preserves the monotonicity of the control-point values, thus avoiding the spurious oscillations
that a cubic b spline can introduce.

In the prior config file, these free parameters are name-keyed as :code:`chem_<species>_mr_<i>`
(species symbol in lower case), running from :code:`chem_<species>_mr_0` at the bottom of the
atmosphere to :code:`chem_<species>_mr_<k-1>` at its top.

Details on the implementation of this parametrisation are discussed
:ref:`here <sec:profile_parametrisations>`.


.. _sec:chemistry_model_sf:

Step function
-------------

The step-function model describes the vertical profile of a single species by two constant
mixing ratios separated by a sharp transition in pressure: a constant deep mixing ratio below
the transition pressure and a constant upper mixing ratio above it. This can, for example,
approximate a species that is quenched in the deep atmosphere but depleted or enhanced higher
up. It is chosen in the ``forward_model.toml`` file with the keyword :code:`sf`
(or :code:`step_function`), followed by the species' formula. This model takes exactly one
parameter:

.. code:: toml

   chemistry = [
     { model = "sf", params = ["CH4"] },
   ]

As for the other free-profile models, only a single species can be used here. The step-function
model always has three free parameters (with the species symbol written in lower case):

  - :code:`chem_<species>_mr_deep`: the constant mixing ratio at pressures greater than or equal
    to the transition pressure, i.e. in the deep atmosphere;

  - :code:`chem_<species>_mr_upper`: the constant mixing ratio at pressures lower than the
    transition pressure, i.e. in the upper atmosphere;

  - :code:`chem_<species>_p_transition`: the transition pressure that separates the two levels.

The transition is a sharp step, without any smoothing: at every level BeAR assigns the deep
mixing ratio if the pressure is at or above the transition pressure, and the upper mixing ratio
otherwise.

As with all chemistry models, these parameters are name-keyed and can be listed in any order in
the prior config file.


.. _sec:chemistry_model_eq:

Equilibrium chemistry
---------------------

BeAR has the option of calculating the chemical composition of the atmosphere using the equilibrium
chemistry code FastChem. In the ``forward_model.toml`` file, this is chosen by

.. code:: toml

   chemistry = [
     { model = "eq", params = ["fastchem_parameters.dat"] },
   ]

The indicated file after the keyword :code:`eq` is the main configuration file for FastChem 4.
It contains information on the elemental abundances and thermochemical data that FastChem
should use. This file has the following structure

.. include:: ../examples/fastchem_parameters_fc4.dat
   :literal:

The parameter file consists of six blocks, in this order:

  1. the elemental abundance file;

  2. the gas-phase species data file. Optionally, a condensate species data file can be given
     on the *same* line, separated by a space (e.g. :code:`logK.dat logK_condensates.dat`), to
     enable equilibrium condensation;

  3. the relative accuracy of the chemistry iteration;

  4. the accuracy of the element conservation (new in FastChem 4);

  5. the maximum number of chemistry iterations;

  6. the maximum number of internal solver iterations.

BeAR already includes samples of the elemental abundance and species data files in the folder
``fastchem_data``. Otherwise, they can be obtained from the
`FastChem repository <https://github.com/NewStrangeWorlds/FastChem>`_. More details on FastChem
and these parameters can be found in the
`FastChem documentation <https://newstrangeworlds.github.io/FastChem/>`_.

The equilibrium chemistry model always has one free parameter, the metallicity factor
:code:`fc_metallicity`. It is a general factor that multiplies all elemental abundances read from the
abundance file, *except* hydrogen and helium. It refers to the elemental abundances specified in the
FastChem parameter file, which are not necessarily always solar.

In addition, individual elemental abundance ratios can be turned into free parameters by listing them
as ``"X/Y"`` strings after the parameter file in the :code:`params` array:

.. code:: toml

   chemistry = [
     { model = "eq", params = ["fastchem_parameters.dat", "C/O", "Fe/H", "V/H"] },
   ]

Each ratio adds a free parameter with the name :code:`fc_<x>_<y>_ratio` (both element symbols in lower
case), e.g. :code:`fc_c_o_ratio` for C/O or :code:`fc_v_h_ratio` for V/H. For a ratio relative to
hydrogen (X/H), the parameter scales the abundance relative to the reference X/H value; for a ratio
X/Y with Y other than hydrogen, the parameter is the literal number ratio of the two elements. As with
all priors, these parameters are name-keyed and can be listed in any order in the prior config file,
with :code:`fc_metallicity` as the metallicity factor.

Even though FastChem is a very fast chemistry code, the computational time will increase substantially when using this chemistry model.
One way to decrease the calculation time is to reduce the number of elements treated in FastChem to only the most important ones, such
as H, He, C, O, and N. Information on how to change the FastChem input files can again be found in the
`FastChem documentation <https://newstrangeworlds.github.io/FastChem/>`_.


.. _sec:chemistry_model_bg:

Background species
------------------

This chemistry model will fill up the background of the atmosphere with a specific chemical species.
The mixing ratio :math:`x_\mathrm{bg}` of this background species :math:`i` is determined by the
mixing ratios of all other species :math:`j`:

.. math::
  x_\mathrm{bg} = 1 - \sum_{j, j \neq i} x_j  \ .

For a gas giant or a brown dwarf, the background is usually given by a combination of molecular
hydrogen (H2) and helium (He), whereas for the atmosphere of an Earth-like planet, the
background species is rather molecular nitrogen (N2).

This background chemistry model is chosen in the ``forward_model.toml`` file by using
the keyword :code:`bg` or :code:`background`, followed by a single parameter that selects one of
three modes:

.. code:: toml

   chemistry = [
     { model = "bg", params = ["N2"] },
   ]

The three available modes are:

  - **A single chemical species symbol** (e.g. :code:`"N2"`): the background is filled up with that
    one species. The species has to be present in the header file ``src/chemistry/chem_species.h``.

  - **"H2He"**: the background is filled with a mixture of H2 and He at their fixed solar elemental
    abundance ratio. This special "species" does not need to be added to the list of chemical species
    in ``src/chemistry/chem_species.h``, but H2 and He individually need to be present there.

  - **"HHeEquilibrium"**: the background is filled with hydrogen and helium, where the split between
    atomic H, molecular H2, and He is computed from a temperature- and pressure-dependent equilibrium.
    This is useful for hot atmospheres in which hydrogen is partly dissociated.

Since the mixing ratio of the background species is determined by those of all other species,
none of these modes requires a free parameter in the prior distribution file.

.. _sec:chemistry_model_mixing:

Mixing different chemistry models
---------------------------------

BeAR also has the ability to use multiple chemistry models simultaneously. For example, to
constrain the constant abundances of some species and then fill up the rest of the atmosphere
with H2 and He, the following can be used:

.. code:: toml

   chemistry = [
     { model = "iso", params = ["H2O", "NH3", "CO2"] },
     { model = "bg", params = ["H2He"] },
   ]

Another use case could be to perform a retrieval with constant mixing ratios for a selection
of species and a separate species, say CH4, that is assumed to have non-constant abundances.
Together with a background gas, the following can be used as configuration
in the ``forward_model.toml`` file:

.. code:: toml

   chemistry = [
     { model = "iso", params = ["H2O", "NH3", "CO2"] },
     { model = "free_cs", params = ["CH4", "5"] },
     { model = "background", params = ["N2"] },
   ]

BeAR will call the chemistry models in the order they appear in this list. That means, in this
example it would first set the constant mixing ratios of H2O, NH3, and CO2, then use a variable
profile based on cubic splines for methane, and finally fill up the rest of the atmosphere with
molecular nitrogen.

The mean molecular weight is recalculated after each chemistry model. BeAR will also check the
sum of all mixing ratios and reject parameter combinations where the sum exceeds 1.1.

Yet another possibility is to calculate the background atmosphere in chemical equilibrium and try
to retrieve a separate species assumed not to be in equilibrium:

.. code:: toml

   chemistry = [
     { model = "eq", params = ["fastchem_parameters.dat"] },
     { model = "iso", params = ["CO"] },
   ]

Here, BeAR would first use FastChem to determine the chemical composition in equilibrium and then replace the CO abundance
by a constant mixing ratio that is a separate free parameter. This can, for example, simulate the impact of vertical mixing.
An extension of that case could also include other non-constant species to take into account photochemical effects in the
upper atmosphere:

.. code:: toml

   chemistry = [
     { model = "eq", params = ["fastchem_parameters.dat"] },
     { model = "iso", params = ["CO"] },
     { model = "free_cs", params = ["SO2", "5"] },
     { model = "free_cs", params = ["HCN", "5"] },
   ]

Because the priors are name-keyed, the free parameters can be listed in any order in the prior
config file, independently of the order of the chemistry models in ``forward_model.toml``.

Thus, for the first example, the following priors need to be listed:

  - constant mixing ratios :code:`chem_h2o_mr`, :code:`chem_nh3_mr`, :code:`chem_co2_mr`

  - 5 mixing ratios :code:`chem_ch4_mr_0` to :code:`chem_ch4_mr_4` for CH4

For the last example above, the prior file needs to list

  - :code:`fc_metallicity` for the equilibrium chemistry

  - the constant CO mixing ratio :code:`chem_co_mr`

  - 5 mixing ratios :code:`chem_so2_mr_0` to :code:`chem_so2_mr_4` for SO2

  - 5 mixing ratios :code:`chem_hcn_mr_0` to :code:`chem_hcn_mr_4` for HCN


One should still pay careful attention to the order of the chemistry models themselves, since they
are applied in sequence. For example, in this case:

.. code:: toml

   chemistry = [
     { model = "iso", params = ["CO"] },
     { model = "eq", params = ["fastchem_parameters.dat"] },
   ]

CO would first be set to an isoprofile and then in the second step be replaced with the results from the equilibrium chemistry
calculation. The free mixing ratio of CO will, therefore, remain unconstrained.

Another incorrect order would, for example, also be the following case:

.. code:: toml

   chemistry = [
     { model = "bg", params = ["H2He"] },
     { model = "iso", params = ["H2O", "CO2", "NH3"] },
   ]

Here, BeAR would first fill up the entire atmosphere with H2 and He. This means that when H2O, CO2,
and NH3 are added through the second model, the sum of all mixing ratios would exceed 1.1.
All parameter combinations during the retrieval calculations would, therefore, be rejected.


.. _sec:chemical_species:

Chemical species list
---------------------

For internal and external communications in the code, BeAR defines a list of all chemical
species it can use. The corresponding source code can be found in the
header file ``src/chemistry/chem_species.h``.

Here, an enumeration is first declared

..  code-block:: cpp
    :caption: src/chemistry/chem_species.h

    enum chemical_species_id {
      _TOTAL, _H, _He, _C, _O, _Fe, _Fep, _Ca, _Ti, _Tip, _H2, _H2O, _CO2, _CO, _CH4, _HCN, _NH3, _C2H2,
      _N2, _Na, _K, _H2S, _Hm, _TiO, _VO, _FeH, _SH, _MgO, _AlO, _CaO, _CrH, _MgH, _CaH, _TiH, _OH, _e,
      _V, _Vp, _Mn, _Si, _Cr, _Crp, _SiO, _SiO2, _SO2, _CS2, _Co, _Ni, _Mg, _13C16O, _SiS, _SiH, _CH, _HCl, _ScH};

This creates a list of integer constants, such that later in the code, the position of, say, H2O
in a vector with chemical species can be referred to as :code:`_H2O`. As can be seen above, the
current list contains about 55 species, including neutral metals (Fe, Ca, Ti, V, Mn, Si, Cr, Co,
Ni, Mg), their ions (Fe+, Ti+, V+, Cr+), free electrons (e-), a wide range of molecules
(CO, CH4, HCN, NH3, TiO, VO, FeH, SO2, SiS, HCl, ...), and the isotopologue 13C16O.

Each species is described by the struct :code:`chemistry_data`, which has four fields: the integer
:code:`id`, the usual :code:`symbol`/formula, the :code:`fastchem_symbol` used by FastChem, and the
atomic/molecular weight:

..  code-block:: cpp
    :caption: src/chemistry/chem_species.h

    struct chemistry_data{
      chemical_species_id id;
      std::string symbol;
      std::string fastchem_symbol;
      double molecular_weight;
    };

A vector with the important chemical information of all species is then created:

..  code-block:: cpp
    :caption: src/chemistry/chem_species.h

    const std::vector<chemistry_data> species_data{
      {_TOTAL,  "Total",  "Total",       0.0},
      {_H,      "H",      "H",           1.00784},
      {_He,     "He",     "He",          4.002602},
      {_C,      "C",      "C",           12.0107},
      {_O,      "O",      "O",           15.999},
      {_Fe,     "Fe",     "Fe",          55.845},
      {_Fep,    "Fe+",    "Fe1+",        55.845},
      {_Ca,     "Ca",     "Ca",          40.078},
      {_Ti,     "Ti",     "Ti",          47.867},
      {_Tip,    "Ti+",    "Ti1+",        47.867},
      {_H2,     "H2",     "H2",          2.01588},
      {_H2O,    "H2O",    "H2O1",        18.01528},
      {_CO2,    "CO2",    "C1O2",        44.01},
      {_CO,     "CO",     "C1O1",        28.0101},
      {_CH4,    "CH4",    "C1H4",        16.04246},
      {_HCN,    "HCN",    "C1H1N1_hcn",  27.0253},
      {_NH3,    "NH3",    "H3N1",        17.03052},
      {_C2H2,   "C2H2",   "C2H2",        26.04},
      {_N2,     "N2",     "N2",          28.0134},
      {_Na,     "Na",     "Na",          22.98977},
      {_K,      "K",      "K",           39.0983},
      {_H2S,    "H2S",    "H2S1",        34.09099},
      {_Hm,     "H-",     "H1-",         1.00784},
      {_TiO,    "TiO",    "O1Ti1",       63.8664},
      {_VO,     "VO",     "O1V1",        66.9409},
      {_FeH,    "FeH",    "H1Fe1",       56.853},
      {_SH,     "SH",     "H1S1",        34.08},
      {_MgO,    "MgO",    "Mg1O1",       40.3044},
      {_AlO,    "AlO",    "Al1O1",       42.981},
      {_CaO,    "CaO",    "Ca1O1",       56.0774},
      {_CrH,    "CrH",    "Cr1H1",       54.0040},
      {_MgH,    "MgH",    "H1Mg1",       26.3209},
      {_CaH,    "CaH",    "Ca1H1",       41.0859},
      {_TiH,    "TiH",    "H1Ti1",       48.87484},
      {_OH,     "OH",     "H1O1",        17.008},
      {_e,      "e-",     "e-",          5.4857990907e-4},
      {_V,      "V",      "V",           50.9415},
      {_Vp,     "V+",     "V1+",         50.9415},
      {_Mn,     "Mn",     "Mn",          54.938044},
      {_Si,     "Si",     "Si",          28.085},
      {_Cr,     "Cr",     "Cr",          51.996},
      {_Crp,    "Cr+",    "Cr1+",        51.996},
      {_SiO,    "SiO",    "O1Si1",       44.08},
      {_SiO2,   "SiO2",   "O2Si1",       60.08},
      {_SO2,    "SO2",    "O2S1",        64.066},
      {_CS2,    "CS2",    "C1S2",        76.139},
      {_Co,     "Co",     "Co",          58.9332},
      {_Ni,     "Ni",     "Ni",          58.6934},
      {_Mg,     "Mg",     "Mg",          24.3050},
      {_13C16O, "13C16O", "C1O1",        28.9982},
      {_SiS,    "SiS",    "S1Si1",       60.150},
      {_SiH,    "SiH",    "H1Si1",       29.093},
      {_CH,     "CH",     "C1H1",        13.01},
      {_HCl,    "HCl",    "Cl1H1",       36.46},
      {_ScH,    "ScH",    "H1Sc1",       47.98}};

This vector connects the defined integer constants with actual chemical species. The
second column contains the usual symbol/formula of the chemical species, the third
the corresponding Hill notation for FastChem, and the last column the atomic/molecular
weight. Note that some FastChem symbols carry a suffix to disambiguate isomers, such as
HCN whose FastChem symbol is :code:`C1H1N1_hcn`.
The molecular-weight information is used internally to compute important quantities, such as the
mean molecular weight. A special species is :code:`_TOTAL`. This refers to the total
number density, given through the pressure via the ideal gas law. It has a molecular
weight of 0 to prevent it from contributing to the mean molecular weight.

BeAR already includes a list of important species. If a species that should be used
in a retrieval is missing, it needs to be added to ``src/chemistry/chem_species.h``.

This has to be done at two places:

  - the enumeration :code:`enum chemical_species_id`

  - the vector :code:`const std::vector<chemistry_data> species_data`

It is recommended to always add a new species at the end of each vector to avoid confusion.
For example, to add phosphine (PH3) as a new species, the code needs to be adapted to

..  code-block:: cpp
    :caption: src/chemistry/chem_species.h

    enum chemical_species_id {
      _TOTAL, _H, _He, _C, _O, /* ... */, _HCl, _ScH, _PH3};

    const std::vector<chemistry_data> species_data{
      {_TOTAL,  "Total",  "Total",       0.0},
      /* ... */
      {_ScH,    "ScH",    "H1Sc1",       47.98},
      {_PH3,    "PH3",    "H3P1",        33.99758}};

After changing the header file, the code needs to be re-compiled.
