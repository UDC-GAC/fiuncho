#!/usr/bin/env bash
set -xe

# Load env vars
source ./.github/vars/$(echo ${GH_MATRIX_CC} | cut -d'+' -f1).sh
source ./.github/vars/$(echo ${GH_MATRIX_CC} | cut -d'+' -f2).sh

# Define CMake arguments
CMAKE_ARGS="-DCMAKE_BUILD_TYPE=RelWithDebInfo \
-DGTABLE_OP_WIDTH=${GH_MATRIX_OP_WIDTH} \
-DMI_OP_WIDTH=${GH_MATRIX_OP_WIDTH}"
SRC_DIR=$(pwd)

# Build
mkdir build/
cd build/
cmake ${CMAKE_ARGS} ${SRC_DIR}
make all
