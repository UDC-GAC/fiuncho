===================================
Technical documentation
===================================

-----------------------------------
Using Fiuncho as a library
-----------------------------------

If you are a developer and want to use any of the C++ classes and methods
implemented in this program, you can do so by linking to the main library
``libfiuncho`` declared in the CMake project.

The following example shows a short program that makes use of Fiuncho's
methods, and its corresponding CMake project configuration. Assume the following
project structure::

   myproject/
   ├── build/
   ├── fiuncho/
   ├── CMakeLists.txt
   └── main.cpp

The directory ``fiuncho/`` contains the source code of Fiuncho and ``build/`` is
the directory where we will build the project from its sources. The program
``main.cpp`` simply reads a data set using the :cpp:class:`Dataset` class
provided in Fiuncho, and prints the number of SNPs read:

.. code-block:: cpp
   :linenos:

   #include <fiuncho/dataset/Dataset.h>
   #include <iostream>

   int main(int argc, char **argv)
   {
      // Arguments to the program are the tped and tfam input files
      if (argc < 3){
         return 1;
      }
      const std::string tped(argv[1]);
      const std::string tfam(argv[2]);
      // Read data
      const auto dataset = Dataset<uint64_t>::read(tped, tfam);
      // Print the number of SNPs read
      std::cout << dataset.snps << std::endl;
      // Exit
      return 0;
   }

The ``CMakeLists.txt`` file will declare the project, include fiuncho as a CMake
subdirectory, declare the main executable target and link it to ``libfiuncho``:

.. code-block:: cmake
   :linenos:

   cmake_minimum_required(VERSION 3.11)
   project(Myproject VERSION 1 LANGUAGES CXX)
   set(CMAKE_CXX_STANDARD 14)
   # Add fiuncho project as a subdirectory
   add_subdirectory(fiuncho)
   # Add our main executable
   add_executable(main main.cpp)
   # Link it against fiuncho's library
   target_link_libraries(main libfiuncho)

``libfiuncho``'s build can be configured using the same CMake variables defined
in the :ref:`installation-and-usage:Advanced configuration` subsection.

-----------------------------------
Classes documentation
-----------------------------------

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Dataset class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: Dataset
   :members:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
GenotypeTable class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: GenotypeTable
   :members:

Function :cpp:func:`GenotypeTable::combine_and_popcnt` has multiple
implementations:

* File ``src/avx512vpopcntdq/gt_popcnt.cpp``:

  .. doxygenfile:: src/avx512vpopcntdq/gt_popcnt.cpp
     :sections: func user-defined

* File ``src/avx512bw/gt_popcnt_avx512bw_hs.cpp``:

  .. doxygenfile:: src/avx512bw/gt_popcnt_avx512bw_hs.cpp
     :sections: func

* File ``src/avx512bw/gt_popcnt_avx512bw_lu.cpp``:

  .. doxygenfile:: src/avx512bw/gt_popcnt_avx512bw_lu.cpp
     :sections: func

* File ``src/avx512bw/gt_popcnt_avx2_cpu.cpp``:

  .. doxygenfile:: src/avx512bw/gt_popcnt_avx2_cpu.cpp
     :sections: func

* File ``src/avx512bw/gt_popcnt_avx2_hs.cpp``:

  .. doxygenfile:: src/avx512bw/gt_popcnt_avx2_hs.cpp
     :sections: func

* File ``src/avx512bw/gt_popcnt_avx2_lu.cpp``:

  .. doxygenfile:: src/avx512bw/gt_popcnt_avx2_lu.cpp
     :sections: func

* File ``src/avx512bw/gt_popcnt_avx2_lu_orig.cpp``:

  .. doxygenfile:: src/avx512bw/gt_popcnt_avx2_lu_orig.cpp
     :sections: func

* File ``src/avx512bw/gt_popcnt_native_movdq.cpp``:

  .. doxygenfile:: src/avx512bw/gt_popcnt_native_movdq.cpp
     :sections: func

* File ``src/avx512bw/gt_popcnt_native_unrolled_errata.cpp``:

  .. doxygenfile:: src/avx512bw/gt_popcnt_native_unrolled_errata.cpp
     :sections: func

* File ``src/avx2/gt_popcnt_avx2_cpu.cpp``:

  .. doxygenfile:: src/avx2/gt_popcnt_avx2_cpu.cpp
     :sections: func

* File ``src/avx2/gt_popcnt_avx2_hs.cpp``:

  .. doxygenfile:: src/avx2/gt_popcnt_avx2_hs.cpp
     :sections: func

* File ``src/avx2/gt_popcnt_avx2_lu.cpp``:

  .. doxygenfile:: src/avx2/gt_popcnt_avx2_lu.cpp
     :sections: func

* File ``src/avx2/gt_popcnt_avx2_lu_orig.cpp``:

  .. doxygenfile:: src/avx2/gt_popcnt_avx2_lu_orig.cpp
     :sections: func

* File ``src/avx2/gt_popcnt_native_movdq.cpp``:

  .. doxygenfile:: src/avx2/gt_popcnt_native_movdq.cpp
     :sections: func

* File ``src/avx2/gt_popcnt_native_unrolled_errata.cpp``:

  .. doxygenfile:: src/avx2/gt_popcnt_native_unrolled_errata.cpp
     :sections: func

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
ContingencyTable class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: ContingencyTable
   :members:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
MutualInformation class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: MutualInformation
   :members:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Distribution class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: Distribution
   :members:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
MPIEngine class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: MPIEngine
   :members:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Search class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: Search
   :members:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
ThreadedSearch class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: ThreadedSearch
   :members:
