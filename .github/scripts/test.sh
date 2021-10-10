#!/usr/bin/env bash
set -xe

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
sde64 ${ARCH_FLAG} -- make test
