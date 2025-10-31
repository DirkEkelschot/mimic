#!/bin/bash

#========================Build parmetis======================================

ROOT=${PWD}

MPICC_PATH=$(which mpicc)
MPI_PREFIX=${MPICC_PATH%/bin/mpicc}
echo "MPI build that is used: $MPI_PREFIX"

#========================Build metis-5.1.0======================================

curl -O https://karypis.github.io/glaros/files/sw/metis/metis-5.1.0.tar.gz
tar -xvf metis-5.1.0.tar.gz
cd metis-5.1.0
mkdir build
cd build
cmake -DCMAKE_VERBOSE_MAKEFILE=1 \
      -DGKLIB_PATH=${ROOT}/metis-5.1.0/GKlib \
      -DCMAKE_INSTALL_PREFIX=${ROOT}/metis-5.1.0/build \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5 ..
make
make install
cd ../../

if [ -f "machine.cmake" ]; then
    echo "set(DEFAULT_METIS-5.1.0_ROOT ${ROOT}/metis-5.1.0/build)" >> machine.cmake
else
    echo "set(DEFAULT_MPI_ROOT ${MPI_PREFIX})" >> machine.cmake
    echo "set(DEFAULT_METIS-5.1.0_ROOT ${ROOT}/metis-5.1.0/build)" >> machine.cmake
fi

