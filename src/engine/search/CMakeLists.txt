# MPI3SNP ~ https://github.com/chponte/mpi3snp
#
# Copyright 2018 Christian Ponte
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
# documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
# Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

add_library(engine_search Search.cpp Search.h)
target_include_directories(engine_search INTERFACE .)
target_link_libraries(engine_search PUBLIC engine_interface utils PRIVATE MPI::MPI_CXX dataset definitions utils)
if (${TARGET_ARCH} STREQUAL "GPU")
    target_link_libraries(engine_search PRIVATE engine_gpu)
else ()
    target_link_libraries(engine_search PRIVATE engine_cpu)
endif ()
