Installation
============

Obtaining the code
~~~~~~~~~~~~~~~~~~

BeAR is hosted on the NewStrangeWorlds GitHub Teams page:
https://github.com/newstrangeworlds/bear. If ``git`` is available,
the repository can be simply cloned with

.. code:: bash

   git clone https://github.com/newstrangeworlds/bear

Prerequisites
~~~~~~~~~~~~~

BeAR is written in C++ and CUDA. It uses features of the C++11
standard and, therefore, requires a compiler that implements this
standard. Furthermore, it needs the CUDA compiler and framework
from NVIDIA.
Additional libraries and codes
(FastChem, MultiNest, CDisort, Boost)
will be downloaded from their respective repositories and be compiled as well.

The complete list of prerequisites for a basic ``CMake`` installation
is:

-  a C++ compiler (e.g. ``g++``)

-  NVIDIA's CUDA framework and drivers

- ``CMake``, at least version 3.10

-  an OpenMP library

The C++ compiler will be detected by the CMake script when it
generates the makefiles. The codes and libraries downloaded and compiled during
the installation of BeAR have addiational requirements.
This includes, in particular,

- a FORTRAN compiler (e.g. ``gfortran``) for the MultiNest library

- a C compiler (e.g. ``gcc``) for the CDisort code


.. _sec:install_config:

Configuration and compilation with CMake
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Before BeAR can be compiled, ``CMake`` is required to
  configure the compilation files, locate libraries, and write the
  makefiles that will perform the actual compilations. If required
  libraries are missing, ``CMake`` will report a corresponding error
  message. In this case, the missing libraries or compilers need to be
  installed before the configuration can be completed.
| To run the ``CMake`` configuration, first create the ``build`` folder
  inside the BeAR source code folder and switch to the folder:

.. code:: bash

   mkdir build
   cd build

For a basic installation, within the folder run ``CMake``:

.. code:: bash

   cmake ..

After ``CMake`` successfully configured the compilation files,
BeAR can be compiled by running:

.. code:: bash

   make

Upon successful compilation, the executable ``bear`` should be
present in the main BeAR folder.
