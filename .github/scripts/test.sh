#!/usr/bin/env bash
set -xe

# Intel MPI needs to load vars for mpirun to be on PATH
source ./.github/vars/$(echo ${GH_MATRIX_CC} | cut -d'+' -f1).sh

# Run all tests using Intel SDE emulating the icelake-server processor
if [ "${GH_MATRIX_TARGET_ARCH}" == "icelake-server" ]; then
    ARCH_FLAG="-icx"
elif [ "${GH_MATRIX_TARGET_ARCH}" == "cascadelake" ]; then
    ARCH_FLAG="-clx"
elif [ "${GH_MATRIX_TARGET_ARCH}" == "haswell" ]; then
    ARCH_FLAG="-hsw"
else
    echo "Target arch not contemplated in the test script"
    exit 1
fi
cd build/
sde64 ${ARCH_FLAG} -- env CTEST_OUTPUT_ON_FAILURE=1 make test
