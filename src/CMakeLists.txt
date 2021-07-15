# This file is part of Fiuncho.
#
# Fiuncho is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Fiuncho is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Fiuncho. If not, see <http://www.gnu.org/licenses/>.

################################# Dependencies #################################

find_package(Threads)
if (NOT TARGET Threads::Threads)
    message(FATAL_ERROR "Missing threading library")
endif()

find_package(MPI)
if (NOT TARGET MPI::MPI_CXX)
    message(FATAL_ERROR "Missing MPI library")
endif()

find_package(AVX)
find_package(VectorLib)

################################# Definitions  #################################

file(GLOB_RECURSE HEADER_LIST "${PROJECT_SOURCE_DIR}/include/fiuncho/*.h")
set(SOURCE_LIST_BASE
    "${PROJECT_SOURCE_DIR}/src/cpu/base/GenotypeTable.cpp"
    "${PROJECT_SOURCE_DIR}/src/cpu/base/ContingencyTable.cpp"
    "${PROJECT_SOURCE_DIR}/src/cpu/base/MutualInformation.cpp")
set(SOURCE_LIST_AVX2
    "${PROJECT_SOURCE_DIR}/src/cpu/avx2/GenotypeTable.cpp"
    "${PROJECT_SOURCE_DIR}/src/cpu/avx2/ContingencyTable.cpp"
    "${PROJECT_SOURCE_DIR}/src/cpu/avx2/MutualInformation.cpp")
set(SOURCE_LIST_AVX512F256
    "${PROJECT_SOURCE_DIR}/src/cpu/avx2/GenotypeTable.cpp"
    "${PROJECT_SOURCE_DIR}/src/cpu/avx2/ContingencyTable.cpp"
    "${PROJECT_SOURCE_DIR}/src/cpu/avx512f256/MutualInformation.cpp")
set(SOURCE_LIST_AVX512F512
    "${PROJECT_SOURCE_DIR}/src/cpu/avx512f512/GenotypeTable.cpp"
    "${PROJECT_SOURCE_DIR}/src/cpu/avx512f512/ContingencyTable.cpp"
    "${PROJECT_SOURCE_DIR}/src/cpu/avx512f512/MutualInformation.cpp")

################################### Targets  ###################################

if(FORCE_AVX512F512)
    if(NOT AVX512F_ENABLED)
        message(FATAL_ERROR "FORCE_AVX512F512 is ${FORCE_AVX512F512}, "
            "but AVX512F is not available")
    endif()
    if (NOT (SVML_AVAILABLE OR LIBMVEC_AVAILABLE))
        message(FATAL_ERROR "FORCE_AVX512F512 is ${FORCE_AVX512F512}, "
            "but there is no math vector library available")
    endif()
    add_library(libfiuncho ${SOURCE_LIST_AVX512F512})
    target_compile_options(libfiuncho PUBLIC "-DALIGN=64")
elseif(FORCE_AVX512F256)
    if(NOT AVX512F_ENABLED)
        message(FATAL_ERROR "FORCE_AVX512F256 is ${FORCE_AVX512F256}, "
            "but AVX512F is not available")
    endif()
    if (NOT (SVML_AVAILABLE OR LIBMVEC_AVAILABLE))
        message(FATAL_ERROR "FORCE_AVX512F256 is ${FORCE_AVX512F256}, "
            "but there is no math vector library available")
    endif()
    add_library(libfiuncho ${SOURCE_LIST_AVX512F256})
    target_compile_options(libfiuncho PUBLIC "-DALIGN=32")
elseif(FORCE_AVX2)
    if(NOT AVX2_ENABLED)
        message(FATAL_ERROR "FORCE_AVX2 is ${$FORCE_AVX2}, "
            "but AVX2 is not available")
    endif()
    if (NOT (SVML_AVAILABLE OR LIBMVEC_AVAILABLE))
        message(FATAL_ERROR "FORCE_AVX2 is ${FORCE_AVX2}, "
            "but there is no math vector library available")
    endif()
    add_library(libfiuncho ${SOURCE_LIST_AVX2})
    target_compile_options(libfiuncho PUBLIC "-DALIGN=32")
elseif(FORCE_NOAVX)
    add_library(libfiuncho ${SOURCE_LIST_BASE})
    target_compile_options(libfiuncho PUBLIC)
else()
    # Default behaviour: use AVX2 if its available
    if (AVX2_ENABLED AND (SVML_AVAILABLE OR LIBMVEC_AVAILABLE))
        add_library(libfiuncho ${SOURCE_LIST_AVX2})
        target_compile_options(libfiuncho PUBLIC "-DALIGN=32")
    else()
        add_library(libfiuncho ${SOURCE_LIST_BASE})
        target_compile_options(libfiuncho PUBLIC)
    endif()
endif()

target_include_directories(libfiuncho PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(libfiuncho PUBLIC Threads::Threads MPI::MPI_CXX)
