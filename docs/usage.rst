==========================================
Usage
==========================================

Fiuncho can only be used in a command-line environment. The program input and
output are provided through files.


------------------------------------------
Command-line usage
------------------------------------------

Fiuncho can be invoked as follows::

   fiuncho [-h] [--version] [-n <integer>]
           [-t <integer>] -o <integer>
           files...


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

files...
    **Required.** List of strings indicating the path of every input and output
    file. Input files go first in any order, output file goes last.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following command executes fiuncho using two processes with 16 threads each,
running a fourth-order analysis. The program reads the input from ``data.tped``
and ``data.tfam``, and writes the top 100 combinations to the file
``output.txt``:

.. code-block:: bash

    mpiexec -n 2 --bind-to numa fiuncho -t 16 -o 4 \
        -n 100 data.tped data.tfam output.txt


------------------------------------------
Data format
------------------------------------------

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Input
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fiuncho can read the input data in two different formats: ``tped`` files
(accompanied by the ``tfam`` file) and ``raw`` files. It uses the file extension
to determine the format used.

TPED file format
""""""""""""""""

Fiuncho can read PLINK ``tped`` and ``tfam`` file formats, representing variants
and their genotype calls, and the different samples' information, respectively.
The complete specification for these formats is available at `PLINK's
documentation <https://www.cog-genomics.org/plink/1.9/formats#tped>`__. Fiuncho
only uses the genotype calls from the ``tped`` file and the phenotype value of
each sample from the ``tfam`` file. The rest of the information is ignored,
although it must be present in the input data.

* **TPED file**

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

* **TFAM file**

  A ``tfam`` file indicates for each sample its family id, within-family id,
  within-family id of father, within-family id of mother, sex and phenotype
  value. The following ``tfam`` file example includes this information for the
  previous eight samples show in the ``tped`` example:

  .. code-block:: plain

      case0 case0 0 0 0 2
      case1 case1 0 0 0 2
      case2 case2 0 0 0 2
      case3 case3 0 0 0 2
      control0 control0 0 0 0 1
      control1 control1 0 0 0 1
      control2 control2 0 0 0 1
      control3 control3 0 0 0 1

RAW file format
""""""""""""""""

The ``raw`` file format represents, in a single file, all variant information
for every sample, as well as the case or control class for every smaple. The
first row of the the file contains the header line, naming the different column.
The following rows contain a sample per row, starting with the sample
information, and followed by the genotype calls (encoded as the minor allele
count for each locus) and phenotype class. The complete specification for this
format is available at `PLINK's documentation
<https://www.cog-genomics.org/plink/1.9/formats#raw>`__.

Fiuncho relies on the header file to identify which columns to read.
Conventional ``raw`` files start with five columns with the sample information
(named ``FID``, ``IID``, ``PAT``, ``MAT`` and ``SEX``) followed by the phenotype
column (named ``PHENOTYPE``, with value ``1`` for controls and ``2`` for cases).
Then, a variable number of columns follow, with two columns per variant (named
``<Variant ID>_{A,C,G,T}`` and ``<Variant ID>_HET``), or one column per variant
if the dominant component (column ``<Variant ID>_HET``) is ommited. Fiuncho will
ignore the dominant component information, regardless if its present or not. The
following ``raw`` file example shows four variants for eight samples:

.. code-block:: plain

    FID      IID      PAT MAT SEX PHENOTYPE N0_A N1_A N2_A N3_A
    case0    case0    0   0   0   2         1    0    0    0
    case1    case1    0   0   0   2         1    0    0    1
    case2    case2    0   0   0   2         0    0    0    0
    case3    case3    0   0   0   2         2    0    1    0
    control0 control0 0   0   0   1         0    0    1    1
    control1 control1 0   0   0   1         2    1    0    0
    control2 control2 0   0   0   1         0    0    0    1
    control3 control3 0   0   0   1         0    0    0    0


Some less conventional ``raw`` files, such as those generated by simulators like
`GAMETES <https://sourceforge.net/projects/gametes/>`_, provide the simulated
data in the ``raw`` format with three main differences:

1. Column names are different.
2. Sample information is absent.
3. The phenotype column (named ``Class``) uses ``0`` to represent controls and
   ``1`` for cases.

Fiuncho will take notice of this and read the information acordingly.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fiuncho provides a list of variant combinations and their associated Mutual
Information as the output. The following file shows a sample output when
searching for fourth-order interactions:

.. code-block:: plain

    0 1 2 11 0.000319958
    0 1 8 11 0.000310183
    0 1 4 11 0.000308275
    0 1 9 11 0.000300407
    0 1 6 11 0.000268698
    0 1 5 11 0.000248909
    0 1 3 11 0.000245333
    0 1 7 11 0.000204086
    2 3 4 5 0.000158548
    4 6 7 8 0.00015223