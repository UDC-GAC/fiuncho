export MPI_PACKAGES="intel-oneapi-mpi intel-oneapi-mpi-devel"
export I_MPI_CC=${CC}
export I_MPI_CXX=${CXX}
export MPICC=mpiicc
export MPICXX=mpiicpc

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
