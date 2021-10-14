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
if (NOT Threads_FOUND)
    message(FATAL_ERROR "Missing threading library")
endif()

find_package(MPI)
if (NOT MPI_FOUND)
    message(FATAL_ERROR "Missing MPI library")
endif()

find_package(AVX)
find_package(VectorLib)

############################### Helper functions  ##############################

function(PREP_GTABLE_512)
    # Check if AVX512BW is supported
    if (NOT AVX512BW_ENABLED)
        message(FATAL_ERROR "Specified GT_OP_WIDTH of 512, but the "
            "AVX512BW vector extension is not available")
    endif()
    # Add genotype table constructor and combination function
    set(SOURCE_LIST "${PROJECT_SOURCE_DIR}/src/avx512bw/gt_combine.cpp")
    # Select contingency table calculation implementation
    if(POPCNT_IMPL)
        if(POPCNT_IMPL STREQUAL "popcnt-512")
            # Check if AVX512_VPOPCNTDQ is supported
            if (NOT AVX512VPOPCNTDQ_ENABLED)
                message(FATAL_ERROR "Specified POPCNT_IMPL=popcnt-512, but
                    the AVX512VPOPCNTDQ vector extension is not available")
            endif()
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx512vpopcntdq/gt_popcnt.cpp")
        elseif(POPCNT_IMPL STREQUAL "harley-seal-512")
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx512bw/gt_popcnt_avx512bw_hs.cpp")
        elseif(POPCNT_IMPL STREQUAL "lookup-512")
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx512bw/gt_popcnt_avx512bw_lu.cpp")
        elseif(POPCNT_IMPL STREQUAL "cpu-256")
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx512bw/gt_popcnt_avx2_cpu.cpp")
        elseif(POPCNT_IMPL STREQUAL "harley-seal-256")
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx512bw/gt_popcnt_avx2_hs.cpp")
        elseif(POPCNT_IMPL STREQUAL "lookup-original-256")
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx512bw/gt_popcnt_avx2_lu_orig.cpp")
        elseif(POPCNT_IMPL STREQUAL "lookup-256")
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx512bw/gt_popcnt_avx2_lu.cpp")
        elseif(POPCNT_IMPL STREQUAL "popcnt-movdq-64")
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx512bw/gt_popcnt_native_movdq.cpp")
        elseif(POPCNT_IMPL STREQUAL "popcnt-unrolled-errata-64")
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx512bw/gt_popcnt_native_unrolled_errata.cpp")
        else()
            message(FATAL_ERROR "Incorrect POPCNT_IMPL function: "
                "${POPCNT_IMPL}")
        endif()
    else()
        if(AVX512VPOPCNTDQ_ENABLED)
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx512vpopcntdq/gt_popcnt.cpp")
        else()
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx512bw/gt_popcnt_avx512bw_lu.cpp")
        endif()
    endif()
    # Update SOURCE_LIST variable in parent scope
    set(SOURCE_LIST ${SOURCE_LIST} PARENT_SCOPE)
    # Define data alignment
    set(DATA_ALIGN 64 PARENT_SCOPE)
endfunction()

function(PREP_GTABLE_256)
    # Check if AVX2 is supported
    if (NOT AVX2_ENABLED)
        message(FATAL_ERROR "Specified GT_OP_WIDTH of 256, but the "
            "AVX2 vector extension is not available")
    endif()
    # Add genotype table constructor and combination function
    set(SOURCE_LIST "${PROJECT_SOURCE_DIR}/src/avx2/gt_combine.cpp")
    # Select contingency table calculation implementation
    if(POPCNT_IMPL)
        if(POPCNT_IMPL STREQUAL "cpu-256")
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx2/gt_popcnt_avx2_cpu.cpp")
        elseif(POPCNT_IMPL STREQUAL "harley-seal-256")
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx2/gt_popcnt_avx2_hs.cpp")
        elseif(POPCNT_IMPL STREQUAL "lookup-256")
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx2/gt_popcnt_avx2_lu.cpp")
        elseif(POPCNT_IMPL STREQUAL "lookup-original-256")
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx2/gt_popcnt_avx2_lu_orig.cpp")
        elseif(POPCNT_IMPL STREQUAL "popcnt-movdq-64")
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx2/gt_popcnt_native_movdq.cpp")
        elseif(POPCNT_IMPL STREQUAL "popcnt-unrolled-errata-64")
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx2/gt_popcnt_native_unrolled_errata.cpp")
        else()
            message(FATAL_ERROR "Incorrect POPCNT_IMPL function: "
                "${POPCNT_IMPL}")
        endif()
    else()
        # By default use the avx lookup algorithm
        list(APPEND SOURCE_LIST
            "${PROJECT_SOURCE_DIR}/src/avx2/gt_popcnt_avx2_lu.cpp")
    endif()
    # Update SOURCE_LIST variable in parent scope
    set(SOURCE_LIST ${SOURCE_LIST} PARENT_SCOPE)
    # Define data alignment
    set(DATA_ALIGN 32 PARENT_SCOPE)
endfunction()

function(PREP_GTABLE_64)
    # Add genotype table constructor and combination function
    set(SOURCE_LIST
        "${PROJECT_SOURCE_DIR}/src/base/gt_combine.cpp"
        "${PROJECT_SOURCE_DIR}/src/base/gt_popcnt.cpp"
        PARENT_SCOPE)
endfunction()

function(PREP_MI_512)
    # Check if AVX512BW is supported
    if (NOT AVX512BW_ENABLED)
        message(FATAL_ERROR "Specified GT_OP_WIDTH of 512, but the "
            "AVX512BW vector extension is not available")
    endif()
    # Check if a math vector library is available
    if (NOT (LIBMVEC_AVAILABLE OR SVML_AVAILABLE))
        message(FATAL_ERROR "Specified GT_OP_WIDTH of 512, but there is no "
            "vector math library")
    endif()
    # Add source files
    list(APPEND SOURCE_LIST
        "${PROJECT_SOURCE_DIR}/src/avx512bw/mi_compute_512.cpp")
    # Define data alignment, overwrite if necessary
    set(DATA_ALIGN 64 PARENT_SCOPE)
    # Update SOURCE_LIST variable in parent scope
    set(SOURCE_LIST ${SOURCE_LIST} PARENT_SCOPE)
endfunction()

function(PREP_MI_256)
    # Check if a math vector library is available
    if (NOT (LIBMVEC_AVAILABLE OR SVML_AVAILABLE))
        message(FATAL_ERROR "Specified GT_OP_WIDTH of 256, but there is no "
            "vector math library")
    endif()
    # Select MI implementation
    if(MI_IMPL)
        if(MI_IMPL STREQUAL "if-nomask")
            # Check if AVX2 is supported
            if (NOT AVX2_ENABLED)
                message(FATAL_ERROR "Specified GT_OP_WIDTH of 256, but the "
                    "AVX2 vector extension is not available")
            endif()
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx2/mi_compute.cpp")
        elseif(MI_IMPL STREQUAL "if-mask")
            # Check if AVX512BW is supported
            if (NOT AVX512BW_ENABLED)
                message(FATAL_ERROR "Specified GT_OP_WIDTH of 512, but the "
                    "AVX512BW vector extension is not available")
            endif()
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx512bw/mi_compute_256.cpp")
        else()
            message(FATAL_ERROR "Incorrect MI_IMPL function: ${MI_IMPL}")
        endif()
    else()
        # If AVX512BW is available, use the AVX512BW implementation of the
        # 256-bit compute operation. Otherwise, use the AVX2 implementation
        if(AVX512BW_ENABLED)
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx512bw/mi_compute_256.cpp")
        else()
            list(APPEND SOURCE_LIST
                "${PROJECT_SOURCE_DIR}/src/avx2/mi_compute.cpp")
        endif()
    endif()
    # Update SOURCE_LIST variable in parent scope
    set(SOURCE_LIST ${SOURCE_LIST} PARENT_SCOPE)
    # Define data alignment only if its not defined (can only be 64, which
    # should be kept; 32, which is the same as here, or undefined)
    if(NOT DATA_ALIGN)
        set(DATA_ALIGN 32 PARENT_SCOPE)
    endif()
endfunction()

function(PREP_MI_64)
    # Add contingency table and MI constructors, and MI compute function
    list(APPEND SOURCE_LIST "${PROJECT_SOURCE_DIR}/src/base/mi_compute.cpp")
    # Update SOURCE_LIST variable in parent scope
    set(SOURCE_LIST ${SOURCE_LIST} PARENT_SCOPE)
endfunction()

################################# Definitions  #################################

file(GLOB_RECURSE HEADER_LIST "${PROJECT_SOURCE_DIR}/include/fiuncho/*.h")

if(GT_OP_WIDTH)
    if(GT_OP_WIDTH EQUAL 512)
        PREP_GTABLE_512()
    elseif(GT_OP_WIDTH EQUAL 256)
        PREP_GTABLE_256()
    elseif(GT_OP_WIDTH EQUAL 64)
        PREP_GTABLE_64()
    else()
        message(FATAL_ERROR "Incorrect GT_OP_WIDTH of ${GT_OP_WIDTH}")
    endif()
else()
    # By default use the highest width available
    if(AVX512BW_ENABLED)
        PREP_GTABLE_512()
    elseif(AVX2_ENABLED)
        PREP_GTABLE_256()
    else()
        PREP_GTABLE_64()
    endif()
endif()

if(MI_OP_WIDTH)
    if(MI_OP_WIDTH EQUAL 512)
        PREP_MI_512()
    elseif(MI_OP_WIDTH EQUAL 256)
        PREP_MI_256()
    elseif(MI_OP_WIDTH EQUAL 64)
        PREP_MI_64()
    else()
        message(FATAL_ERROR "Incorrect MI_OP_WIDTH of ${GT_OP_WIDTH}")
    endif()
else()
    # By default use the highest width available
    if (LIBMVEC_AVAILABLE OR SVML_AVAILABLE)
        if(AVX512BW_ENABLED)
            PREP_MI_512()
        elseif(AVX2_ENABLED)
            PREP_MI_256()
        else()
            PREP_MI_64()
        endif()
    else()
        PREP_MI_64()
    endif()
endif()

################################### Targets  ###################################

add_library(libfiuncho ${SOURCE_LIST})
if(DATA_ALIGN)
    target_compile_options(libfiuncho PUBLIC "-DALIGN=${DATA_ALIGN}")
endif()
target_include_directories(libfiuncho PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(libfiuncho PUBLIC Threads::Threads MPI::MPI_CXX)
