**Install MIMIC**
  1. Get the repository:
      git clone https://developer.nasa.gov/dekelsch/mimic.git

  2. Modify machine.cmake and set appropriate paths to the local install directories for MPI, Metis, ParMetis, Boost and TinyXml.
     Additionally, paths can be set to HDF5_DIR and VTK_DIR in order to find preinstalled versions of these dependencies. 
	
     cp -r machine_NAS.cmake machine.cmake when installing MIMIC on NASA Advanced Supercomputing (NAS) facility.

  3. By default MIMIC downloads and installs automatically Mmg and ParMmg HDF5 and VTK (as CMake external projects):
      cd mimic
      mkdir build
      cd build
      cmake ..
      make
      make install
