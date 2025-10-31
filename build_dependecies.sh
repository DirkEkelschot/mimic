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
ln -s ${ROOT}/gklib/gklib-install/lib64/libGKlib.a ${ROOT}/gklib/gklib-install/lib/libGKlib.a
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

#========================Build boost c++ libraries==============================

curl -O https://archives.boost.io/release/1.83.0/source/boost_1_83_0.tar.gz
tar -xvf boost_1_83_0.tar.gz
cd boost_1_83_0
mkdir boost-install
./bootstrap.sh --prefix=${ROOT}/boost_1_83_0/boost-install
./b2
./b2 install

touch machine.cmake
echo "set(DEFAULT_MPI_ROOT ${MPI_PREFIX})" > machine.cmake
echo "set(DEFAULT_METIS_ROOT ${ROOT}/metis/metis-install)" >> machine.cmake
echo "set(DEFAULT_PARMETIS_ROOT ${ROOT}/parmetis/build)" >> machine.cmake
echo "set(DEFAULT_METIS-5.1.0_ROOT ${ROOT}/metis-5.1.0/build)" >> machine.cmake
echo "set(DEFAULT_BOOST_ROOT ${ROOT}/boost_1_83_0/boost-install)" >> machine.cmake
