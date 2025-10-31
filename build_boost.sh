#!/bin/bash

#========================Build parmetis======================================

ROOT=${PWD}

MPICC_PATH=$(which mpicc)
MPI_PREFIX=${MPICC_PATH%/bin/mpicc}
echo "MPI build that is used: $MPI_PREFIX"


curl -O https://archives.boost.io/release/1.83.0/source/boost_1_83_0.tar.gz
tar -xvf boost_1_83_0.tar.gz
cd boost_1_83_0
mkdir boost-install
./bootstrap.sh --prefix=${ROOT}/boost_1_83_0/boost-install
./b2
./b2 install

MACHINE_FILE="${ROOT}/machine.cmake"

if [ -f "$MACHINE_FILE" ]; then
    echo "set(DEFAULT_BOOST_ROOT ${ROOT}/boost_1_83_0/boost-install)" >> "$MACHINE_FILE"
else
    echo "set(DEFAULT_MPI_ROOT ${MPI_PREFIX})" >> "$MACHINE_FILE"
    echo "set(DEFAULT_BOOST_ROOT ${ROOT}/boost_1_83_0/boost-install)" >> "$MACHINE_FILE"
fi

