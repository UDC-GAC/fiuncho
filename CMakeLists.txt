############################### CMake variables  ###############################

cmake_minimum_required(VERSION 3.11)
cmake_policy(SET CMP0074 NEW)
project(Fiuncho VERSION 1 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)
set(CMAKE_STATIC_LIBRARY_PREFIX "")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

############################# Build configuration  #############################

# Compiler configurations:
# GNU
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG
        "-DDEBUG -g -fverbose-asm -Wall -Wextra -Wpedantic")
    set(CMAKE_CXX_FLAGS_RELEASE
        "-DNDEBUG")
    set(CMAKE_CXX_FLAGS_BENCHMARK
        "-DBENCHMARK -g")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
        "-g")
# Intel
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS_DEBUG
        "-DDEBUG -g -debug inline-debug-info -Wall -qopt-report=5")
    set(CMAKE_CXX_FLAGS_RELEASE
        "-DNDEBUG")
    set(CMAKE_CXX_FLAGS_BENCHMARK
        "-DBENCHMARK -g")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
        "-g")
# Clang
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS_DEBUG
        "-DDEBUG -g -fverbose-asm -Wall -Wextra -Wpedantic")
    set(CMAKE_CXX_FLAGS_RELEASE
        "-DNDEBUG")
    set(CMAKE_CXX_FLAGS_BENCHMARK
        "-DBENCHMARK -g")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
        "-g")
else()
    message(FATAL_ERROR "Unsupported compiler detected: "
        ${CMAKE_CXX_COMPILER})
endif()

# Build type configuration
set(DEFAULT_BUILD_TYPE "RELEASE")
set_property(GLOBAL PROPERTY DEBUG_CONFIGURATIONS "DEBUG;RELWITHDEBINFO")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS
        "Setting build type to ${DEFAULT_BUILD_TYPE} as none was specified")
    set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}"
        CACHE STRING "Choose the type of build" FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE
        PROPERTY STRINGS "DEBUG" "RELEASE" "BENCHMARK" "RELWITHDEBINFO")
else()
    string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE)
endif ()

############################# Project directories  #############################

# Library
add_subdirectory(src)
# If this is the top-level project, add more targets
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    # Project documentation
    add_subdirectory(docs)
    # Main program
    add_subdirectory(app)
    # Tests
    enable_testing()
    add_subdirectory(tests)
    # Benchmarks
    add_subdirectory(benchmarks EXCLUDE_FROM_ALL)
endif()
