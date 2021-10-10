#!/usr/bin/env bash
set -xe

# Load env vars
source ./.github/vars/$(echo ${GH_MATRIX_CC} | cut -d'+' -f1).sh
source ./.github/vars/$(echo ${GH_MATRIX_CC} | cut -d'+' -f2).sh

# Install deps
sudo apt-get update -qq 1>/dev/null
sudo apt-get upgrade -qq 1>/dev/null
sudo apt-get install ${COMPILER_PACKAGES} ${MPI_PACKAGES} -qq 1>/dev/null

# Install Intel SDE
INTEL_SDE_URL="https://software.intel.com/content/dam/develop/external/\
us/en/documents/downloads/sde-external-8.69.1-2021-07-18-lin.tar.bz2"
INTEL_SDE_TARFILE=$(echo ${INTEL_SDE_URL} | rev | cut -d '/' -f 1 | rev)
wget -q ${INTEL_SDE_URL}
mkdir intel_sde/
tar -xf ${INTEL_SDE_TARFILE} --strip-components=1 -C intel_sde/
echo "$(pwd)/intel_sde" >> $GITHUB_PATH
