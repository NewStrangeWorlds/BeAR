
.. _sec:chemistry_models:

Chemistry Models
=================

BeAR currently includes the following chemistry models:

  - :ref:`Equilibrium chemistry <sec:chemistry_model_eq>`

  - :ref:`Background species <sec:chemistry_model_bg>`

  - :ref:`Isoprofiles <sec:chemistry_model_iso>`
  
  - :ref:`Isoprofiles with centred-log-ratio priors <sec:chemistry_model_iso_clr>`

  - :ref:`Piecewise polynomials <sec:chemistry_model_poly>`

  - :ref:`Cubic b splines <sec:chemistry_model_cs>`

Details on the implementation of the latter four parametrisations are 
discussed :ref:`here <sec:profile_parametrisations>`. 

BeAR can also mix several different of these models in a single retrieval. This
is explained :ref:`in this section <sec:chemistry_model_mixing>`. All chemical 
species that BeAR uses have to be defined in a corresponding header file.
This is described at the :ref:`bottom <sec:chemical_species>` of this chapter. 


.. _sec:chemistry_model_iso:

Isoprofiles
-----------

In the corresponding ``forward_model.config`` of the chosen forward model, 
the isoprofile chemistry model is chosen by using the keyword :code:`iso` 
(or :code:`isoprofile`),  followed by a list of species that should be used by the retrieval:

.. code:: 

   #Retrieved chemical species
   iso  H2O CH4 NH3 K H2S CO2

For each species, a corresponding entry in the prior config file is required.
The order of the priors has to be identical to the one in the ``forward_model.config``.

Furthermore, as described by `Kitzmann et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020ApJ...890..174K/>`_,
the abundance of sodium (Na) is linked to the one of potassium (K) through the ratio of their
solar elemental abundances. Thus, in the example above, the Na abundance would automatically be
derived from K, without explicitly including the former as a free parameter.

If Na should be retrieved independently, it needs to be added as a separate species as shown below.

.. code:: 

   #Retrieved chemical species
   iso  H2O CH4 NH3 K H2S CO2 Na
   

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
produce :math:`\xi_j` values subject to the constraints that :math:`\min\left(\mathbf{x}\right) = 10^{-10}` and :math:`\max\left(\mathbf{x}\right) = 1`. 
It is important to note that the prior boundaries for :math:`\xi_j` depend on the number of molecules in the retrieval and the 
chosen value of the smallest allowed mixing ratio. For uniform priors between 0 and 1, BeAR will do the conversion to :math:`\xi_j` and then
to the mixing ratios :math:`x_j` automatically.

This chemistry model is chosen in the ``forward_model.config`` file by using

.. code:: 

   #Retrieved chemical species
   iso_clr H2O CH4 NH3 H2S CO2 

As noted above, for these 5 species only 4 priors are needed due to the constraint that the sum of their mixing ratios has to be unity.
This chemistry model will, therefore, select one of these free species as the dominant background gas. In the prior config file
only priors for the first four molecules would need to be listed. The abundance of the last species in the list above
will be determined from the other molecules.
   

.. _sec:chemistry_model_poly:   

Piecewise polynomials
---------------------

If a chemical species, say, H2O should be included using non-constant abundances by using
the parametrisation with piecewise polynomials, the keyword :code:`free`, followed by the species'
formula, the number of elements :code:`k`, and the polynomial degree :code:`q` need to be added to 
the ``forward_model.config`` file:

.. code:: 

   #Retrieved chemical species
   free H2O k q

Note, that unlike the isoprofile chemistry above, only a single species can be used here.
As stated :ref:`here <sec:profile_parametrisations>`, the number of free parameters for this
parametrisation is :math:`k q + 1`.
The free parameters that have to be added to the priors config file correspond to the mixing ratios 
of the chosen species at the :math:`k q + 1` discrete pressure points that are distributed 
equidistantly in logarithmic pressure space. The first parameter in that list refers to the bottom 
of the atmosphere, while the last represents its top.


.. _sec:chemistry_model_cs:   

Cubic b splines
---------------

For parametrising the vertical profile of a species, for example, CH4 by a cubic b spline, the 
keyword :code:`free_cs` (or :code:`free_cspline`), followed by the species'
formula and the number of points :code:`k` need to be added to  the ``forward_model.config`` file:

.. code:: 

   #Retrieved chemical species
   free_cs CH4 k

Note, that like for the piecewise polynomials, only a single species can be used here.
As stated :ref:`here <sec:profile_parametrisations>`, the number of free parameters for this
parametrisation is equal to the number of points :math:`k` and needs to be at least 5.
They correspond to the mixing ratios at the :math:`k` discrete pressure points that are distributed
equidistantly in logarithmic pressure space.

These free parameters have to be added to the priors config file in descending pressures. 
The first parameter in the prior list refers to the bottom of the atmosphere, while the last represents its top.


.. _sec:chemistry_model_eq:

Equilibrium chemistry
---------------------

BeAR has the option of calculating the chemical composition of the atmosphere using the equilibrium
chemistry code FastChem. In the ``forward_model.config`` file, this is chosen by 

.. code:: 

   #Retrieved chemical species
   eq fastchem_parameters.dat
   
The indicated file after the keyword :code:`eq` is the main configuration file for FastChem.
It contains information on the elemental abundances and thermochemical data that FastChem
should use. This file has the following structure

.. include:: ../examples/fastchem_parameters.dat
   :literal:
   
The first and second entry correspond to the location of the elemantal abundance file and the file with
the equilibrium constants, respectively. BeAR already includes a sample of these files in the folder
``fastchem_data``. Otherwise, they can be obtained from the `FastChem repository <https://github.com/NewStrangeWorlds/FastChem>`_.

The third entry is the relative accuracy of the FastChem convergence criterion, followed by the limits for the number of chemistry and 
solver iterations. More details on FastChem and the parameters can be found in the 
`FastChem documentation <https://newstrangeworlds.github.io/FastChem/>`_.

The equilibrium chemistry model has two free parameters: the metallicity factor :code:`M/H` and the :code:`C/O` ratio. :code:`M/H` is
a general factor that the read-in elemental abundances are multiplied with. Thus, both :code:`M/H` and the :code:`C/O` parameters
refer to the elemental abundances that are specified in the FastChem parameter file, which are not necessarily always solar. Both have to
be added to the prior configuration file, with :code:`M/H` as the first and the :code:`C/O` ratio as the second parameter.

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

This background chemistry model is chosen in the ``forward_model.config`` file by using
the keyword :code:`bg` or :code:`background`

.. code:: 

   #Retrieved chemical species
   bg N2

followed by the chemical species that will fill up the atmosphere. For the special case
of a background mixture of both H2 and He with their solar elemental abundance ratio, the special
species code :code:`H2He` is used. This special "species" does not need to be added to the
list of chemical species in the header file ``src/chemistry/chem_species.h`` but rather
H2 and He need to be present there individually.

Since the mixing ratio of the background species is determined by those of all other species,
no free parameter is required for this model in the prior distribution file.

.. _sec:chemistry_model_mixing:

Mixing different chemistry models
---------------------------------

BeAR also has the ability to use multiple chemistry models simultaneously. For example, to
constrain the constant abundances of some species and then fill up the rest of the atmosphere
with H2 and He, the following can be used:

.. code:: 

   #Retrieved chemical species
   iso H2O NH3 CO2
   bg H2He

Another use case could be to perform a retrieval with constant mixing ratios for a selection 
of species and a separate species, say CH4, that is assumed to have non-constant abundances.
Together with a background gas, the following can be used as configuration 
in the ``forward_model.config`` file:

.. code:: 

   #Retrieved chemical species
   iso H2O NH3 CO2
   free_cs CH4 5
   background N2

BeAR will call the chemistry models in the order they appear in this list. That means, in this 
example it would first set the constant mixing ratios of H2O, NH3, and CO2, then use a variable 
profile based on cubic splines for methane, and finally fill up the rest of the atmospere with 
molecular nitrogen.

The mean molecular weight is recalculated after each chemistry model. BeAR will also check the 
sum of all mixing ratios and neglect models where the sum exceeds unity.

Yet another possibility is to calculate the background atmosphere in chemical equilibrium and try 
to retrieve a separate species assumed not to be in equilibrium:

.. code:: 

   #Retrieved chemical species
   eq fastchem_parameters.dat
   iso CO
   
Here, BeAR would first use FastChem to determine the chemical composition in equilibrium and then replace the CO abundance
by a constant mixing ratio that is a separate free parameter. This can, for example, simulate the impact of vertical mixing.
An extension of that case could also include other non-constant species to take into account photochemistry effects in the
upper atmosphere:

.. code:: 

   #Retrieved chemical species
   eq fastchem_parameters.dat
   iso CO
   free_cs SO2 5
   free_cs HCN 5

In the prior configuration file, the free parameters have to appear in the same order as the chemistry models listed in 
``forward_model.config``. 

Thus, for the first example, the following priors need to be listed:
  
  - constant mixing ratios for H2O, NH3, CO2
  
  - 5 mixing ratios for CH4

For the last example above, the prior file needs to list

  - :code:`M/H` and :code:`C/O` for the equilibrium chemistry
  
  - the constant CO mixing ratio
  
  - 5 mixing ratios for SO2
  
  - 5 mixing ratios for HCN


One should pay careful attention to the order of the chemistry models that should be used. For example, in this case:

.. code:: 

   #Retrieved chemical species
   iso CO
   eq fastchem_parameters.dat
   
CO would first be set to an isoprofile and then in the second step be replaced with the results from the equilibrium chemistry
calculation. The free mixing ratio of CO will, therefore, remain unconstrained.

Another incorrect order would, for example, also be the following case:

.. code:: 

   #Retrieved chemical species
   bg H2He
   iso H2O CO2 NH3

Here, BeAR would first fill up the entire atmosphere with H2 and He. This means that when H2O, CO2,
and NH3 are added through the second model, the sum of all mixing ratios would exceed unity.
All parameter combinations during the retrieval calculations would, therefore, be neglected.


.. _sec:chemical_species:

Chemical species list
---------------------

For internal and external communications in the code, BeAR defines a list of all chemical
species it can use. The corresponding source code can be found in the 
header file ``src/chemistry/chem_species.h``.

Here, an enumeration is first declared

..  code-block:: cpp
    :caption: src/chemistry/chem_species.h
    
    enum chemical_species_id {_TOTAL, _H, _He, _C, _O, _H2, _H2O, _CO2};

This creates a list of integer constants, such that later in the code, the position of, say, H2O 
in a vector with chemical species can be referred to as :code:`_H2O`. Next, a vector
with important chemical information is created:

..  code-block:: cpp
    :caption: src/chemistry/chem_species.h
    
    const std::vector<chemistry_data> species_data{ 
      {_TOTAL, "Total", "Total",  0.0},
      {_H,     "H",     "H",      1.00784},
      {_He,    "He",    "He",     4.002602},
      {_C,     "C",     "C",      12.0107},
      {_O,     "O",     "O",      15.999},
      {_H2,    "H2",    "H2",     2.01588},
      {_H2O,   "H2O",   "H2O1",   18.01528},
      {_CO2,   "CO2",   "C1O2",   44.01}};
    
This vector connects the defined integer constants with actual chemical species. The
second column contains the usual symbol/formula of the chemical species, the third
the corresponding Hill notation for FastChem, and the last column the atomic/molecular
weight.
The latter information is used internally to compute important quantities, such as the
mean molecular weight. A special species is :code:`_TOTAL`. This refers to the total
number density, given through the pressure via the ideal gas law. It has a molecular
weight of 0 to prevent it from contributing to the mean molecular weight.

BeAR already includes a list of important species. If a species that should be used
in a retrieval is missing, it needs to be added to ``src/chemistry/chem_species.h``.

This has to be done at two places:
  
  - the enumeration :code:`enum chemical_species_id`
  
  - the vector :code:`const std::vector<chemistry_data> species_data`

It is recommended to always add a new species at the end of each vector to avoid confusion.
For example, to add methane (CH4) as a new species, the code needs to be adapted to

..  code-block:: cpp
    :caption: src/chemistry/chem_species.h
    
    enum chemical_species_id {_TOTAL, _H, _He, _C, _O, _H2, _H2O, _CO2, _CH4};
    
    const std::vector<chemistry_data> species_data{ 
      {_TOTAL, "Total", "Total",  0.0},
      {_H,     "H",     "H",      1.00784},
      {_He,    "He",    "He",     4.002602},
      {_C,     "C",     "C",      12.0107},
      {_O,     "O",     "O",      15.999},
      {_H2,    "H2",    "H2",     2.01588},
      {_H2O,   "H2O",   "H2O1",   18.01528},
      {_CO2,   "CO2",   "C1O2",   44.01},
      {_CH4,   "CH4",   "C1H4",   16.04246}};

After changing the header file, the code needs to be re-compiled.


