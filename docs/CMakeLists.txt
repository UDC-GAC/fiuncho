################################# Dependencies #################################

# System dependencies
find_package(Sphinx)
if (NOT SPHINX_FOUND)
	message(WARNING "Could not find Sphinx, docs target disabled")
    return()
endif()

find_package(Doxygen)
if (NOT DOXYGEN_FOUND)
	message(WARNING "Could not find Doxygen, docs target disabled")
    return()
endif()

################################ Configuration  ################################

# Configure Doxygen configuration file
set(DOXYGEN_INPUT_DIR "../src ../include")
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
# Doxygen won't create this for us
file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})

################################### Targets  ###################################

# Doxygen main output file
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/xml/index.xml)
# Doxygen target
add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    MAIN_DEPENDENCY ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    COMMENT "Generating doxygen xml"
    VERBATIM)
add_custom_target(doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

# Sphinx main output file
set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)
# Sphinx target
add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
    COMMAND
        ${SPHINX_EXECUTABLE} -b html
        # Tell Breathe where to find the Doxygen output
        -Dbreathe_projects.fiuncho=${DOXYGEN_OUTPUT_DIR}/xml
    ${SPHINX_SOURCE} ${SPHINX_OUTPUT_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS
    # Other docs files you want to track should go here (or in some variable)
    ${CMAKE_CURRENT_SOURCE_DIR}/index.rst
    ${DOXYGEN_INDEX_FILE}
    MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
    COMMENT "Generating documentation with Sphinx")
add_custom_target(sphinx ALL DEPENDS ${SPHINX_INDEX_FILE})

# Docs target
add_custom_target(docs DEPENDS doxygen sphinx)
