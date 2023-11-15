PARMETIS_HOME = /Users/dekelsch/Software/parmetis-4.0.3/install
METIS_HOME = /Users/dekelsch/Software/parmetis-4.0.3/metis/install
HDF5_HOME = /Users/dekelsch/Software/hdf5-1.14.2/install
VTK_HOME = /Users/dekelsch/Software/VTK-9.3.0.rc1/build
MPICH_HOME = /Users/dekelsch/Software/mpich-4.1.2/mpich-4.1.2-install
PARMMG_HOME = /Users/dekelsch/Software/pmmg/build
MMG_HOME = /Users/dekelsch/Software/pmmg/build/Mmg-prefix/src/Mmg-build
XML_HOME = /Users/dekelsch/Software/tinyxml
#CGNS_HOME = /Users/dekelsch/Software/CGNS-3.3.0/src/install
BOOST_HOME = /Users/dekelsch/Software/boost_1_71_0
CXXFLAGS += -std=c++11 -I$(VTK_HOME)/include -I$(MMG_HOME)/include -I$(PARMMG_HOME)/include -I$(PARMETIS_HOME)/include -I$(MPICH_HOME)/include -I$(HDF5_HOME)/include -I$(METIS_HOME)/include -I$(XML_HOME) -I$(BOOST_HOME)

LDFLAGS += -L$(VTK_HOME)/lib -L$(MMG_HOME)/lib -L$(PARMMG_HOME)/lib -L$(PARMETIS_HOME)/lib -L$(METIS_HOME)/lib -L$(MPICH_HOME)/lib -L$(HDF5_HOME)/lib -L$(XML_HOME)

LDLIBS += -lmetis -lparmetis -lhdf5 -lmpi -llapack -lblas -lmmg -lparmmg

CC = /Users/dekelsch/Software/mpich-4.1.2/mpich-4.1.2-install/bin/mpic++

