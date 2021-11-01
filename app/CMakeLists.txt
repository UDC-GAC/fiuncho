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
