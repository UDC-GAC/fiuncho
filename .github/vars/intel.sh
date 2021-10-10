export COMPILER_PACKAGES="intel-basekit intel-hpckit cmake"
export CC=icc
export CXX=icpc
export CFLAGS="-march=${GH_MATRIX_TARGET_ARCH} -mtune=${GH_MATRIX_TARGET_ARCH} \
-O3 -fp-model=fast -Wall"
export CXXFLAGS="${CFLAGS}"

# Add Intel repo
APT_FILE=/etc/apt/sources.list.d/oneAPI.list
if [ ! -f ${APT_FILE} ]; then
    wget -q "https://apt.repos.intel.com/intel-gpg-keys/\
GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB"
    sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
    rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
    echo "deb https://apt.repos.intel.com/oneapi all main" | sudo tee \
        ${APT_FILE}
fi
