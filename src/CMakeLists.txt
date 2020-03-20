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

# System dependencies
find_package(Threads)
if (NOT TARGET Threads::Threads)
    return()
endif()

find_package(MPI)
if (NOT TARGET MPI::MPI_CXX)
    return()
endif()

# External dependencies
include(FetchContent)

FetchContent_Declare(args
    GIT_REPOSITORY https://github.com/Taywee/args
    GIT_TAG 6.2.2)
FetchContent_GetProperties(args)
if(NOT args_POPULATED)
    FetchContent_Populate(args)
    add_subdirectory(${args_SOURCE_DIR} ${args_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

################################# Definitions  #################################

file(GLOB_RECURSE HEADER_LIST "${PROJECT_SOURCE_DIR}/include/fiuncho/*.h")
set(SOURCE_LIST 
    "${PROJECT_SOURCE_DIR}/src/dataset/Dataset.cpp"
    "${PROJECT_SOURCE_DIR}/src/engine/Search.cpp"
    "${PROJECT_SOURCE_DIR}/src/engine/cpu/ContingencyTable.cpp"
    "${PROJECT_SOURCE_DIR}/src/engine/cpu/Cpu_node_information.cpp"
    "${PROJECT_SOURCE_DIR}/src/engine/cpu/CPUEngine.cpp"
    "${PROJECT_SOURCE_DIR}/src/engine/cpu/MutualInformation.cpp"
    "${PROJECT_SOURCE_DIR}/src/utils/Arg_parser.cpp"
    "${PROJECT_SOURCE_DIR}/src/utils/IOMpi.cpp"
    "${PROJECT_SOURCE_DIR}/src/utils/Node_information.cpp"
    "${PROJECT_SOURCE_DIR}/src/utils/Statistics.cpp")

################################### Targets  ###################################

configure_file(${PROJECT_SOURCE_DIR}/include/fiuncho/Definitions.h.in 
    ${PROJECT_BINARY_DIR}/include/fiuncho/Definitions.h)

add_library(fiuncho_lib ${SOURCE_LIST})

target_include_directories(fiuncho_lib PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}/include)

target_link_libraries(fiuncho_lib PUBLIC args Threads::Threads MPI::MPI_CXX)
