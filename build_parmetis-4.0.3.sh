#!/bin/bash

#========================Build parmetis======================================

ROOT=${PWD}

MPICC_PATH=$(which mpicc)
MPI_PREFIX=${MPICC_PATH%/bin/mpicc}
echo "MPI build that is used: $MPI_PREFIX"

curl -O https://karypis.github.io/glaros/files/sw/parmetis/parmetis-4.0.3.tar.gz
tar -xvf parmetis-4.0.3.tar.gz
cd parmetis-4.0.3/metis
mkdir metis-install
make config prefix=${ROOT}/parmetis-4.0.3/metis/metis-install
make
make install

cd ${ROOT}/parmetis-4.0.3
mkdir parmetis-install
make config prefix=${ROOT}/parmetis-4.0.3/parmetis-install
make
make install

MACHINE_FILE="${ROOT}/machine.cmake"

if [ -f "$MACHINE_FILE" ]; then
    echo "set(DEFAULT_METIS_ROOT ${ROOT}/parmetis-4.0.3/metis/metis-install)" >> "$MACHINE_FILE"
    echo "set(DEFAULT_PARMETIS_ROOT ${ROOT}/parmetis-4.0.3/parmetis-install)" >> "$MACHINE_FILE"
else
    echo "set(DEFAULT_MPI_ROOT ${MPI_PREFIX})" >> "$MACHINE_FILE"
    echo "set(DEFAULT_METIS_ROOT ${ROOT}/parmetis-4.0.3/metis/metis-install)" >> "$MACHINE_FILE"
    echo "set(DEFAULT_PARMETIS_ROOT ${ROOT}/parmetis-4.0.3/parmetis-install)" >> "$MACHINE_FILE"
fi


