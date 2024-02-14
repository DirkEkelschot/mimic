PARMETIS_HOME = /Users/dekelsch/Software/parmetis-install
METIS_HOME = /Users/dekelsch/Software/metis-install
HDF5_HOME = /Users/dekelsch/Software/hdf5-1.14.2/install
VTK_HOME = /Users/dekelsch/Software/VTK-9.3.0/build
MPICH_HOME = /Users/dekelsch/Software/mpich-4.1.2
PARMMG_HOME = /Users/dekelsch/Software/pmmg/build
MMG_HOME = /Users/dekelsch/Software/pmmg/build/Mmg-prefix/src/Mmg-build
XML_HOME = /Users/dekelsch/Software/tinyxml
#CGNS_HOME = /Users/dekelsch/Software/CGNS-3.3.0/src/install
BOOST_HOME = /Users/dekelsch/Software/boost_1_83_0
TETGEN_HOME = /Users/dekelsch/Software/tetgen1.6.0
CXXFLAGS += -std=c++11 -I$(VTK_HOME)/include/vtk-9.3 -I$(MMG_HOME)/include -I$(PARMMG_HOME)/include -I$(PARMETIS_HOME)/include -I$(MPICH_HOME)/include -I$(HDF5_HOME)/include -I$(METIS_HOME)/include -I$(XML_HOME) -I$(BOOST_HOME)

LDFLAGS += -L$(VTK_HOME)/lib -L$(MMG_HOME)/lib -L$(PARMMG_HOME)/lib -L$(PARMETIS_HOME)/lib -L$(METIS_HOME)/lib -L$(MPICH_HOME)/lib -L$(HDF5_HOME)/lib -L$(XML_HOME)

LDLIBS += -lmetis -lparmetis -lhdf5 -lmpi -llapack -lblas -lmmg -lparmmg
CC = /Users/dekelsch/Software/mpich-4.1.2/mpich-install/bin/mpic++

