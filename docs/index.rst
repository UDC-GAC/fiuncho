================
Fiuncho
================

|cpp-version| |nbsp| |license| |nbsp| |ctest| |nbsp| |docs|

----------------
What is Fiuncho?
----------------

Fiuncho is an epistasis detection program, implementing an exhaustive search of
epistasis interactions of any given size. Fiuncho uses a `SPMD
<https://en.wikipedia.org/wiki/SPMD>`_ (Single Program, Multiple Data) approach
to accelerate the computation, exploiting the multiple cores available to a
processor, or the multiple processors available in a cluster environment.
Furthermore, Fiuncho takes advantage of the `SIMD
<https://en.wikipedia.org/wiki/SIMD>`_ (Single Instruction, Multiple Data)
vector operations available to each CPU core, offering multiple *AVX*
implementations to support a wide variety of *x86_64* processors.

-----------------------
Contents of this manual
-----------------------

.. toctree::
   :maxdepth: 3

   installation
   usage
   technical-doc
   license

------------------------
Reporting problems
------------------------

If you find any problem using Fiuncho, an error in the documentation or
something that can be improved, please open an `Issue
<https://github.com/UDC-GAC/fiuncho/issues/new>`_ on Github describing the
problem.


.. |nbsp| unicode:: 0xA0
.. |cpp-version| image:: https://img.shields.io/badge/C++-17-blue.svg?style=flat&logo=c%2B%2B
.. |license| image:: https://img.shields.io/github/license/UDC-GAC/fiuncho?color=blue
.. |ctest| image:: https://github.com/UDC-GAC/fiuncho/actions/workflows/ctest.yml/badge.svg?branch=dev
.. |docs| image:: https://readthedocs.org/projects/fiuncho/badge/?version=latest