#!/bin/bash

#========================Build parmetis======================================

ROOT=${PWD}

MPICC_PATH=$(which mpicc)
MPI_PREFIX=${MPICC_PATH%/bin/mpicc}
echo "MPI build that is used: $MPI_PREFIX"
git clone https://github.com/KarypisLab/GKlib.git gklib
cd gklib
mkdir gklib-install
make config cc=mpicc prefix=${ROOT}/gklib/gklib-install
make
make install
cp -r ${ROOT}/gklib/gklib-install/lib64/libGKlib.a ${ROOT}/gklib/gklib-install/lib/libGKlib.a
cd ..

git clone https://github.com/KarypisLab/METIS.git metis
cd metis
mkdir metis-install
make config cc=mpicc gklib_path=${ROOT}/gklib/gklib-install prefix=${ROOT}/metis/metis-install
make
make install
cd ..

git clone https://github.com/KarypisLab/ParMETIS.git parmetis
cd parmetis
mkdir build
cd build
cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
      -DCMAKE_C_COMPILER=mpicc \
      -DGKLIB_PATH=${ROOT}/gklib/gklib-install \
      -DMETIS_PATH=${ROOT}/metis/metis-install \
      -DCMAKE_INSTALL_PREFIX=${ROOT}/parmetis/build ..
make
make install
cd ../../

if [ -f "machine.cmake" ]; then
    echo "set(DEFAULT_METIS_ROOT ${ROOT}/metis/metis-install)" >> machine.cmake
    echo "set(DEFAULT_PARMETIS_ROOT ${ROOT}/parmetis/build)" >> machine.cmake
else
    echo "set(DEFAULT_MPI_ROOT ${MPI_PREFIX})" >> machine.cmake
    echo "set(DEFAULT_METIS_ROOT ${ROOT}/metis/metis-install)" >> machine.cmake
    echo "set(DEFAULT_PARMETIS_ROOT ${ROOT}/parmetis/build)" >> machine.cmake
fi

