export COMPILER_PACKAGES="clang-12 cmake"
export CC=clang-12
export CXX=clang++-12
export CFLAGS="-march=${GH_MATRIX_TARGET_ARCH} -mtune=${GH_MATRIX_TARGET_ARCH} \
-O3 -ffast-math -Wall"
export CXXFLAGS="${CFLAGS}"
