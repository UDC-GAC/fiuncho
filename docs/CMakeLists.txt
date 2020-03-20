# This file is part of Fiuncho.
# Copyright (C) 2020 by Christian Ponte
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

find_package(Sphinx)
if (NOT DEFINED SPHINX_EXECUTABLE)
    return()
endif()

find_package(Doxygen)
if (NOT DEFINED DOXYGEN_EXECUTABLE)
    return()
endif()

include(FetchContent)
FetchContent_Declare(breathe
    GIT_REPOSITORY https://github.com/michaeljones/breathe
    GIT_TAG v4.14.1)
FetchContent_GetProperties(breathe)
if(NOT breathe_POPULATED)
    FetchContent_Populate(breathe)
endif()

################################### Targets  ###################################

# Doxygen
set(DOXYGEN_INPUT_DIR 
    "${Fiuncho_SOURCE_DIR}/src ${Fiuncho_SOURCE_DIR}/include")
set(DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/doxygen")
set(DOXYGEN_XML_OUTPUT_DIR "xml")
configure_file(Doxyfile Doxyfile @ONLY)
add_custom_target(doxygen
    ALL
    COMMAND ${DOXYGEN_EXECUTABLE}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Building code documentation with Doxygen"
    SOURCES Doxyfile)

# Sphynx
configure_file(conf.py conf.py @ONLY)
add_custom_target(docs
    ALL
    COMMAND ${SPHINX_EXECUTABLE} -q -b html -c ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/html
    DEPENDS doxygen
    COMMENT "Building HTML documentation with Sphinx")