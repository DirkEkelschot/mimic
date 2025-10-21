 1) gcc/10.3                4) szip/2.1.1        7) us3d/1.2-dev(default)            10) cuda/12.2.2  
 2) comp-intel/2020.4.304   5) hdf5/1.8.18_mpt   8) mpi-hpe/mpt.2.28_25Apr23_rhel87  
 3) mpi-hpe/mpt.2.30        6) parmetis/4.0.3    9) mpi-hpe/mpt

cmake -DCMAKE_CXX_COMPILER=/nasa/hpe/mpt/2.30_rhel810/bin/mpicxx -DKokkos_DIR=/nobackupp19/dekelsch/Software/kokkos-4.5.01/build -DKokkosKernels_DIR=/nobackupp19/dekelsch/Software/kokkos-kernels/build ..
