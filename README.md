# MPI3SNP

<p>
    <a href="https://doi.org/10.1177/1094342019852128" alt="Publication">
    <img src="https://img.shields.io/badge/DOI-10.1177%2F1094342019852128-blue"/></a>
</p>

## What is MPI3SNP?

MPI3SNP is a parallel software tool dedicated to genome-wide association studies, performing a third-order exhaustive
search. It is targeted to cluster architectures, and mitigates the cubic time complexity inherent to third-order 
searches by exploiting the several layers of parallelism present in a supercomputer. CPU and GPU implementations are
offered.

## Building

Support is currently limited to linux distributions only.

### Requirements

* CMake (>3.0 version)
* A C++14 compatible compiler
* MPI library

Optional:
* CUDA 

### Compilation

CMake is the project build manager. CMake should be able to determine installed compilers and libraries. If this is
is not the case, please refer to your CMake version's documentation. By default, CMake will check for a CUDA
installation and set the target architecture accordingly. This behaviour can be manually controlled by setting the
`TARGET_ARCH` CMake variable to `CPU` or `GPU`. 

Building the sources looks like this:

```commandline
cd MPI3SNP/project/path
mkdir build
cd build/
cmake ..
make -j4
```

## Usage

MPI3SNP takes two files as the input, using the [PLINK/TPED](https://www.cog-genomics.org/plink2/formats#tfam) format,
and writes the results to a third file. All file paths are provided to the program as positional arguments as follows:

```commandline
./MPI3SNP <path/to/tped> <path/to/tfam> <path/to/output>
```

Additional configuration options (specific to either CPU or GPU implementation) are available to the user, and can be
consulted using the `-h` flag.

## Sample files

[Sample files](https://github.com/chponte/mpi3snp/wiki/Sample-files) can be found on MPI3SNP's wiki. These are a syntetic dataset used for performance evaluation, which describe the input file format and can be used for verification/evaluation purposes.

## Troubleshooting

Support is currently limited to linux distributions only. If you are having trouble building/using the application, 
please submit a new issue to get help.

## License

This software is licensed under the GPU GPLv3 license. Check the [LICENSE](LICENSE.md) file for details.
