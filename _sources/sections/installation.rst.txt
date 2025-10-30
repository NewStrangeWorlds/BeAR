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
generates the make files. The codes and libraries downloaded and compiled during
the installation of BeAR have additional requirements.
This includes, in particular,

- a FORTRAN compiler (e.g. ``gfortran``) for the MultiNest library

- the BLAS linear algebra library for MultiNest (e.g. ``OpenBLAS``)

- a C compiler (e.g. ``gcc``) for the CDisort code

BeAR also comes with its own Python interface pyBeAR. The interface allows the forward model
and retrieval calculations to be performed in Python. The interface requires

 - Python 3.9 or higher
 - the pyMultiNest package (if the retrieval should be performed with MultiNest under Python)

 Some combinations of the CUDA, g++, and Python compiler versions seem to produce issues during either
 the compilation or when executing the code. Symptoms include, for example, segmentation faults when 
 the BeAR C++ code is called from within Python. This seems to be caused by compatibility issues within 
 the PyBind11 library and are, thus, not directly related to BeAR and also not fixable within BeAR itself.
 The following combinations have been tested and are known to work with the PyBind11 version used in BeAR:

- g++/gcc/gfortran 12.3
- CUDA 12.0
- Python 3.10


.. _sec:install_config:

Configuration and compilation of BeAR with CMake
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Before BeAR can be compiled, ``CMake`` is required to
  configure the compilation files, locate libraries, and write the
  make files that will perform the actual compilations. If required
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


Configuration and compilation of pyBeAR with CMake
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| To compile the Python interface in addition to BeAR's stand-alone version, 
  The ``CMake`` command above needs to be extended by the following option:

.. code:: bash

   cmake -DUSE_PYTHON=ON ..

After ``CMake`` successfully configured the compilation files,
BeAR can be compiled by running:

.. code:: bash

   make

This will compile the Python interface *and* the stand-alone version of BeAR. Since
the compilation will in fact be run twice, it will obviously take longer than
compiling only the stand-alone version.

After successful compilation, the Python package will be located in the ``python/lib`` 
folder.