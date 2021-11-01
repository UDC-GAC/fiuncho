# Fiuncho

![](https://img.shields.io/badge/C++-17-blue.svg?style=flat&logo=c%2B%2B)
![](https://img.shields.io/github/license/UDC-GAC/fiuncho?color=blue)
[![CTest](https://github.com/UDC-GAC/fiuncho/actions/workflows/ctest.yml/badge.svg)](https://github.com/UDC-GAC/fiuncho/actions/workflows/ctest.yml)
[![Documentation Status](https://readthedocs.org/projects/fiuncho/badge/)](https://fiuncho.readthedocs.io/en/latest/)

Fiuncho is an epistasis detection program, implementing an exhaustive search of
epistasis interactions of any given size. Fiuncho uses a
[SPMD](https://en.wikipedia.org/wiki/SPMD) (Single Program, Multiple Data)
approach to accelerate the computation, exploiting the multiple cores available
to a processor, or the multiple processors available in a cluster environment.
Furthermore, Fiuncho takes advantage of the
[SIMD](https://en.wikipedia.org/wiki/SIMD) (Single Instruction, Multiple Data)
vector operations available to each CPU core, offering multiple *AVX*
implementations to support a wide variety of *x86_64* processors.

## Documentation

Fiuncho's documentation is hosted in RTD:
https://fiuncho.readthedocs.io/en/latest/

## Reporting problems

If you find any problem using Fiuncho, an error in the documentation or
something that can be improved, please open an [Issue][1] describing the problem.

[1]: https://github.com/chponte/fiuncho/issues/new
