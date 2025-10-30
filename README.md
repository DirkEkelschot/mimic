**Install MIMIC**
  1. Get the repository:
      git clone https://developer.nasa.gov/dekelsch/mimic.git

  2. MIMIC relies on MPI, Metis, ParMetis, Boost, ParMMG HDF5 and VTK. The user is required to install MPI, Metis, ParMetis and Boost themselves.
     and ParMMG HDF5 and VTK are installed automatically using the cmake build system. If MIMIC is used for commercial purposes, make sure that 
     the appropriate licenses are in place in order to use ParMetis. ParMetis is free to use for academic and research purposes.
     The following two steps are required to carefully consider. However, in order to simplify them for the user, a build_parmetis.sh script
     can be run to install the ParMetis dependency.

  3. The following steps are required to install Metis and ParMetis

        ** git clone https://github.com/KarypisLab/GKlib.git gklib
        ** cd gklib
        ** mkdir gklib-install
        ** make config cc=mpicc prefix=/path/to/gklib/gklib-install
        ** make
        ** make install

        ** git clone https://github.com/KarypisLab/METIS.git metis
        ** cd metis
        ** mkdir metis-install
        ** make config cc=mpicc gklib_path=/path/to/gklib/gklib-install prefix=/path/to/metis/metis-install
        ** make
        ** make install

        ** git clone https://github.com/KarypisLab/ParMETIS.git parmetis
        ** cd parmetis
        ** mkdir build
	  ** cd build
        ** cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 
                 -DCMAKE_C_COMPILER=mpicc 
                 -DCMAKE_POLICY_VERSION_MINIMUM=3.5 
                 -DGKLIB_PATH=/path/to/gklib/gklib-install 
                 -DMETIS_PATH=/path/to/metis/metis-install 
                 -DCMAKE_INSTALL_PREFIX=/path/to/parmetis/build ..
        ** make
        ** make install

  4. ParMMG relies on Metis but a specific version of Metis (v5.1.0). We therefore install this version of Metis also:

        ** curl -O https://karypis.github.io/glaros/files/sw/metis/metis-5.1.0.tar.gz
        ** tar -xvf metis-5.1.0.tar.gz
        ** cd metis-5.1.0
        ** mkdir build
	  ** cd build
        ** cmake -DCMAKE_VERBOSE_MAKEFILE=1 
                 -DGKLIB_PATH=/path/to/metis-5.1.0/GKlib 
                 -DCMAKE_INSTALL_PREFIX=/path/to/metis-5.1.0/build
                 -DCMAKE_C_COMPILER=mpicc 
                 -DCMAKE_POLICY_VERSION_MINIMUM=3.5 ..
        ** make
        ** make install
   
  If the user runs ./build_parmetis.sh, the machine.cmake file is generated automatically.

  5. Modify machine.cmake and set appropriate paths to the local install directories for MPI, Metis, ParMetis, Boost.
     Additionally, paths can be set to HDF5_DIR and VTK_DIR in order to find preinstalled versions of these dependencies. 
	
     cp -r machine_NAS.cmake machine.cmake when installing MIMIC on NASA Advanced Supercomputing (NAS) facility.

  6. By default MIMIC downloads and installs automatically Mmg and ParMmg HDF5 and VTK (as CMake external projects):
      cd mimic
      mkdir build
      cd build
      cmake ..
      make
      make install
