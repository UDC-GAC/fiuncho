#!/usr/bin/env bash
set -e
# Errors will cause the shell script to exit immediately

targets="westmere haswell skylake-avx512"

# Build Fiuncho for each target architecture
for t in $(echo $targets); do
    mkdir -p ${TRAVIS_BUILD_DIR}/build/$t
    cd ${TRAVIS_BUILD_DIR}/build/$t
    CFLAGS="-O3 -march=$t -mtune=$t" CXXFLAGS=$CFLAGS \
            cmake -DCMAKE_BUILD_TYPE=Release \
            ${TRAVIS_BUILD_DIR}
    make all
    # Run tests only if current CPU supports the instruction set
    if [[ $(./test_mi 1>/dev/null 2>&1 || echo $?) != 132 ]]; then
        make test
    fi
done
