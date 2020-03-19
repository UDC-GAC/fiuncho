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

 cmake_minimum_required(VERSION 3.8)
project(MPI3SNP VERSION 1.0)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)

find_package(Threads REQUIRED)
find_package(MPI REQUIRED)
include(CheckLanguage)
check_language(CUDA)

if (NOT DEFINED TARGET_ARCH)
    if (CMAKE_CUDA_COMPILER)
        set(TARGET_ARCH "GPU" CACHE STRING "Select target architecture")
        enable_language(CUDA)
    else ()
        set(TARGET_ARCH "CPU" CACHE STRING "Select target architecture")
    endif ()
    set_property(CACHE TARGET_ARCH PROPERTY STRINGS "CPU" "GPU")
endif()

# Build type configuration
set(DEFAULT_BUILD_TYPE "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to ${DEFAULT_BUILD_TYPE} as none was specified.")
    set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
            STRING "Choose the type of build" FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif ()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3")
if (${TARGET_ARCH} STREQUAL "GPU")
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 11)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_FLAGS "-arch=sm_30")
    set(CMAKE_CUDA_FLAGS_DEBUG "-g -O0")
    set(CMAKE_CUDA_FLAGS_RELEASE "-DNDEBUG -O3")
#    set(CUDA_PROPAGATE_HOST_FLAGS OFF)
endif ()

# Third party libraries (common to GPU and CPU builds)
add_subdirectory(thirdparty)
# Project sources
add_subdirectory(src)
