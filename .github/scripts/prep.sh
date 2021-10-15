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
wget -q --user-agent="Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 \
(KHTML, like Gecko) Chrome/51.0.2704.103 Safari/537.36" ${INTEL_SDE_URL}
mkdir intel_sde/
tar -xf ${INTEL_SDE_TARFILE} --strip-components=1 -C intel_sde/
echo "$(pwd)/intel_sde" >> $GITHUB_PATH

# Post install actions
bash -c "${CC_POST_INSTALL}"
bash -c "${MPI_POST_INSTALL}"
