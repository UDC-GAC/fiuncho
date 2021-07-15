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

include(FetchContent)
FetchContent_Declare(tclap_git
    GIT_REPOSITORY https://github.com/mirror/tclap
    GIT_TAG 1.4.0-rc1)
FetchContent_GetProperties(tclap_git)
if(NOT tclap_git_POPULATED)
    FetchContent_Populate(tclap_git)
    add_subdirectory(
        ${tclap_git_SOURCE_DIR}
        ${tclap_git_BINARY_DIR}
        EXCLUDE_FROM_ALL)
endif()

################################### Targets  ###################################

# Read git hash
execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
# Create string with compile flags
string(STRIP "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}" COMPILER_FLAGS)
# Add all available targets
add_executable(fiuncho main.cpp)
target_link_libraries(fiuncho TCLAP libfiuncho)
target_compile_definitions(fiuncho PRIVATE
    FIUNCHO_VERSION="v${Fiuncho_VERSION}"
    FIUNCHO_COMMIT_HASH="${GIT_HASH}"
    COMPILER_NAME="${CMAKE_CXX_COMPILER_ID}"
    COMPILER_VERSION="${CMAKE_CXX_COMPILER_VERSION}"
    COMPILER_FLAGS="${COMPILER_FLAGS}")
