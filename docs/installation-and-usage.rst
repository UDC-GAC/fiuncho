==========================================
Installation and usage
==========================================

Fiuncho is distributed as a group of source files that the user is required to
compile into a binary before its use. This section describes the process of
downloading the source code, building the application, and using its
command-line interface. The format used for the input data files is also
documented.

For this manual, we assume a standard Linux environment where all dependencies can be downloaded through the distribution's repositories.

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

.. TIP::
    Older versions of the compiler and system libraries have not been tested and
    may work fine, however the performance obtained may not be optimal.

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

More advanced configurations are also possible through CMake's project
variables. In addition to the default CMake variables, this project introduces
the following variables:

CMAKE_BUILD_TYPE
  The default CMake variable to select a build configuration. Accepted values
  are ``Debug``, ``DebWithRelInfo``, ``Release`` and ``Benchmark``.

FORCE_AVX512F512
  Force CMake to build Fiuncho using the AVX Intrinsics implementation using 512
  bit operations from the ``AVX512F`` extension. Accepted values are ``ON`` and
  ``OFF``. This option is incompatible with any other ``FORCE_*`` option.

FORCE_AVX512F256
  Force CMake to build Fiuncho using the AVX Intrinsics implementation using 256
  bit operations from the ``AVX512F`` extension. Accepted values are ``ON`` and
  ``OFF``. This option is incompatible with any other ``FORCE_*`` option.

FORCE_AVX2
  Force CMake to build Fiuncho using the AVX Intrinsics implementation using 256
  bit operations from the ``AVX2`` extension. Accepted values are ``ON`` and
  ``OFF``. This option is incompatible with any other ``FORCE_*`` option.

FORCE_NOAVX
  Force CMake to build Fiuncho not using any of the AVX Intrinsics
  implementations. Accepted values are ``ON`` and ``OFF``. This option is
  incompatible with any other ``FORCE_*`` option.

------------------------------------------
Command-line usage
------------------------------------------

Fiuncho can be invoked as follows::

   fiuncho [-h] [--version] [-n <integer>]
           [-t <integer>] -o <integer>
           tped tfam output


Note that Fiuncho is an MPI program, and as such, it should be called through
``mpiexec`` or any other parallel job launcher such as ``srun`` from SLURM. If
you need help with launching an MPI program, please refer to the MPI or job
scheduling system documentation instead.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Named arguments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-o, --order
    **Required.** Integer equal or greater than 2 specifying the order of the
    epistasis interactions to explore during the search.

-t, --threads
    An integer greater than 0 indicating the number of threads per process to
    use during the search. Note that if you are running an MPI job with multiple
    processes, each process will create the same number of threads. If it's not
    specified, fiuncho will use as many threads as physical cores are available
    to each process.

-n, --noutputs
    An integer greater than 0 indicating the number of combinations to output.
    If it's not specified, it will output 10 combinations.

-h, --help
    Displays usage information and exits.

--version
    Displays version information and exits.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Positional arguments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

tped
    **Required.** First positional argument indicating the path to the tped data
    file.
tfam
    **Required.** Second positional argument indicating the path to the tfam
    data file.
output
    **Required.** Third positional argument indicating the path to the output
    file.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following command executes fiuncho using two processes with 16 threads each,
running a fourth-order analysis on two input files ``data.tped`` and
``data.tfam``, and writing the top 100 combinations to the file ``output.txt``:

.. code-block:: bash

    mpiexec -n 2 --bind-to numa fiuncho -t 16 -o 4 \
        -n 100 data.tped data.tfam output.txt

------------------------------------------
Input data format
------------------------------------------

Fiuncho uses the PLINK ``tped`` and ``tfam`` file formats to represent variants
and their genotype calls, and the different samples' information, respectively.
The complete specification for these formats is available at `PLINK's
documentation <https://www.cog-genomics.org/plink/1.9/formats>`__. Fiuncho only
uses the genotype calls from the ``tped`` file and the phenotype value of each
sample from the ``tfam`` file. The rest of the information is ignored, although
it must be present in the input data.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
tped file format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To briefly summarize it, ``tped`` files contain no header lines, and each line
represents a different variant with four preceding fields: chromosome code,
variant identifier, chromosome position and base-pair coordinate. After these
four fields, the genotype calls for all samples are included. The following
``tped`` file example shows four variants for eight samples:

.. code-block:: plain

    0 N0 0 0 A C C A C C A A C C A A C C C C
    1 N1 0 0 C C C C C C C C C C A C C C C C
    2 N2 0 0 C C C C C C A C C A C C C C C C
    3 N3 0 0 C C A C C C C C A C C C A C C C

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
tfam file format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A ``tfam`` file indicates for each sample its family id, within-family id,
within-family id of father, within-family id of mother, sex and phenotype value.
The following ``tfam`` file example includes this information for the previous
eight samples show in the ``tped`` example:

.. code-block:: plain

    case0 case0 0 0 0 2
    case1 case1 0 0 0 2
    case2 case2 0 0 0 2
    case3 case3 0 0 0 2
    control0 control0 0 0 0 1
    control1 control1 0 0 0 1
    control2 control2 0 0 0 1
    control3 control3 0 0 0 1
