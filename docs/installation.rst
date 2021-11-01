==========================================
Installation
==========================================

Fiuncho is distributed as a group of source files that the user is required to
compile into a binary before its use. This section describes the process of
downloading the source code, building the application, and using its
command-line interface. The format used for the input data files is also
documented.

For this manual, we assume a standard Linux environment where all dependencies
can be downloaded through the distribution's repositories.

------------------------------------------
Downloading the source files
------------------------------------------

All Fiuncho releases can be downloaded from the Github repository
https://github.com/chponte/fiuncho, listed under the Releases tab.

------------------------------------------
Configuring and building the CMake project
------------------------------------------

Fiuncho uses the CMake build system to configure and compile the final
application.

The requirements to build Fiuncho are:

*  Git.
*  CMake :guilabel:`>=3.11`.
*  A C++ compiler. The supported compilers are:

   + GCC :guilabel:`>=8.3.0` with glibc :guilabel:`>=2.22`.
   + Clang :guilabel:`>=10.0.0` with glibc :guilabel:`>=2.22`.
   + Intel C/C++ Compiler :guilabel:`>=19.0`.

*  An MPI library. The supported MPI libraries are:

   + OpenMPI :guilabel:`>=4.0.0`.
   + MPICH :guilabel:`>=3.2`.
   + Intel MPI :guilabel:`>=19.0`.

.. tip::
    Older versions of the compiler and system libraries have not been tested and
    may work fine, however the performance obtained may not be optimal.

.. note::
    Fiuncho has also been built succesfully in Windows using `MSYS2
    <https://www.msys2.org/docs/what-is-msys2/>`_ :guilabel:`3.2.0-15` and
    `MS-MPI
    <https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi>`_
    :guilabel:`10.1.1-5`. This setup is discouraged for unexperienced users.

If the requirements are met, Fiuncho can be configured and built with the
following commands::

    mkdir build
    cd build
    CFLAGS="-O3 -march=native" CXXFLAGS=$CFLAGS cmake ..
    make fiuncho

Compilation flags can be passed to the compiler using the ``CFLAGS`` and
``CXXFLAGS`` environmental variables, as indicated in the example. Specifying an
optimization level 3 (with ``-O3``) and a target compilation architecture (with
``-march``) is highly recommended. If you are building Fiuncho to be run on a
different machine than the current one, replace the ``native`` target
architecture with the appropriate one.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Advanced configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning::
  It is not recommended to change any of these variables. CMake will select the
  best implementation available based on your target CPU.

More advanced configurations are also possible through CMake's project
variables. In addition to the default CMake variables, this project introduces
the following variables:

CMAKE_BUILD_TYPE
  The default CMake variable to select a build configuration. Accepted values
  are ``Debug``, ``DebWithRelInfo``, ``Release`` and ``Benchmark``.

GT_OP_WIDTH
  CMake variable to select the vector width for the operations used during the
  :cpp:class:`GenotypeTable` computation functions
  :cpp:func:`GenotypeTable::combine` and
  :cpp:func:`GenotypeTable::combine_and_popcnt`. Accepted values are: ``512``
  (default if ``AVX512BW`` is available), ``256`` (default if only ``AVX2`` is
  available), ``64`` (default if no AVX extensions are available).

POPCNT_IMPL
  Select the implementation to be used during the :cpp:class:`ContingencyTable`
  computation function :cpp:func:`GenotypeTable::combine_and_popcnt`. Accepted
  values depend on the vector width used:

  * ``GT_OP_WIDTH`` = ``512``:

    * ``popcnt-512`` (default if ``AVX512VPOPCNTDQ`` is available)
    * ``harley-seal-512``
    * ``lookup-512`` (default if only ``AVX512BW`` is available)
    * ``cpu-256``
    * ``harley-seal-256``
    * ``lookup-original-256``
    * ``lookup-256``
    * ``popcnt-movdq-64``
    * ``popcnt-unrolled-errata-64``

  * ``GT_OP_WIDTH`` = ``256``:

    * ``cpu-256``
    * ``harley-seal-256``
    * ``lookup-original-256``
    * ``lookup-256`` (defaulf if only ``AVX2`` is available)
    * ``popcnt-movdq-64``
    * ``popcnt-unrolled-errata-64``

  The documentation for each of these functions is available in the
  :cpp:class:`GenotypeTable` class documentation.

MI_OP_WIDTH
  Select the vector width for the operations used during the
  :cpp:class:`MutualInformation` computation function
  :cpp:func:`MutualInformation::compute`. Accepted values are: ``512`` (default
  if ``AVX512BW`` is available), ``256`` (default if only ``AVX2`` is
  available), ``64`` (default if no AVX extensions are available).

MI_IMPL
  Select the implementation to be used during the :cpp:class:`MutualInformation`
  computation function :cpp:func:`MutualInformation::compute`. Only available
  for ``MI_OP_WIDTH`` = ``256``. Accepted values are:

  * ``if-nomask`` (default if ``AVX512BW`` is available)
  * ``if-mask`` (default if only ``AVX2`` is available)

  The documentation for the two functions is available in the
  :cpp:class:`MutualInformation` class documentation.
