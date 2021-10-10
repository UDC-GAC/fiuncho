export COMPILER_PACKAGES="gcc-10 g++-10 cmake"
export CC=gcc-10
export CXX=g++-10
export CFLAGS="-march=${GH_MATRIX_TARGET_ARCH} -mtune=${GH_MATRIX_TARGET_ARCH} \
-O3 -ffast-math -Wall"
export CXXFLAGS="${CFLAGS}"
