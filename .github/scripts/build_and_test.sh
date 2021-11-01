#!/usr/bin/env bash
set -xe

# Load env vars
source ./.github/vars/$(echo ${GH_MATRIX_CC} | cut -d'+' -f1).sh
source ./.github/vars/$(echo ${GH_MATRIX_CC} | cut -d'+' -f2).sh

# Declare configurations
declare -A cmake_config=(\
    ["icelake-server"]="\
        512,popcnt-512,512,\
        64,,64,"\
    ["cascadelake"]="\
        512,harley-seal-512,,\
        512,lookup-512,,\
        512,cpu-256,,\
        512,harley-seal-256,,\
        512,lookup-original-256,,\
        512,lookup-256,,\
        512,popcnt-movdq-64,,\
        512,popcnt-unrolled-errata-64,,\
        ,,256,if-mask\
        64,,64,"\
    ["haswell"]="\
        256,cpu-256,,\
        256,harley-seal-256,,\
        256,lookup-256,,\
        256,lookup-original-256,,\
        256,popcnt-movdq-64,,\
        256,popcnt-unrolled-errata-64,,\
        64,,64,")
declare -A sde_arch_flag=(\
    ["icelake-server"]="-icx"\
    ["cascadelake"]="-clx"\
    ["haswell"]="-hsw")
src_dir=$(pwd)

# Build loop
mkdir build/
cd build/
for arch in ${!cmake_config[@]}; do
    for config in $(echo ${cmake_config[${arch}]}); do
        # Split config into vars
        gt_width=$(echo ${config} | cut -d',' -f1)
        popcnt_impl=$(echo ${config} | cut -d',' -f2)
        mi_width=$(echo ${config} | cut -d',' -f3)
        mi_impl=$(echo ${config} | cut -d',' -f4)
        # Compose CMake command
        cmake_args="-DCMAKE_BUILD_TYPE=RelWithDebInfo"
        if [ ! -z "$gt_width" ]; then
            cmake_args+=" -DGT_OP_WIDTH=${gt_width}"
        fi
        if [ ! -z "$popcnt_impl" ]; then
            cmake_args+=" -DPOPCNT_IMPL=${popcnt_impl}"
        fi
        if [ ! -z "$mi_width" ]; then
            cmake_args+=" -DMI_OP_WIDTH=${mi_width}"
        fi
        if [ ! -z "$mi_impl" ]; then
            cmake_args+=" -DMI_IMPL=${mi_impl}"
        fi
        # Build
        CFLAGS="${CFLAGS} -march=${arch} -mtune=${arch}" CXXFLAGS="${CFLAGS}" \
            cmake ${cmake_args} ${src_dir}
        make all
        # Test
        sde64 ${sde_arch_flag[${arch}]} -- \
            env CTEST_OUTPUT_ON_FAILURE=1 make test
        # Clean
        rm -rf *
    done
done
