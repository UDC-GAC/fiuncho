export MPI_PACKAGES="openmpi-bin libopenmpi3 libopenmpi-dev"
export OMPI_CC=${CC}
export OMPI_CXX=${CXX}
export MPICC=mpicc
export MPICXX=mpic++

# Fix openmpi complaining about not having enough slots
OPENMPI_HOSTFILE=/etc/openmpi/openmpi-default-hostfile
export MPI_POST_INSTALL="\
if [ -z '$(cat ${OPENMPI_HOSTFILE} | grep -v -E '^#')' ]; then\
    echo \"localhost slots=32\" | sudo tee -a ${OPENMPI_HOSTFILE};\
fi"
