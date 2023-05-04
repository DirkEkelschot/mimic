PARMETIS_HOME = /Users/dekelsch/Software/parmetis-4.0.3/parmetis-install
METIS_HOME = /Users/dekelsch/Software/parmetis-4.0.3/metis/metis-install
HDF5_HOME = /Users/dekelsch/Software/hdf5-1.12.0/hdf5-install-mpich
MPICH_HOME = /Users/dekelsch/Software/mpich-3.3.1/mpich-3.1.1-install
PARMMG_HOME = /Users/dekelsch/Software/parmmg_0805/build
MMG_HOME = /Users/dekelsch/Software/mmg_0805/build
XML_HOME = /Users/dekelsch/Software/tinyxml
CGNS_HOME = /Users/dekelsch/Software/CGNS-4.0.0/build/src
BOOST_HOME = /Users/dekelsch/Software/boost_1_71_0

CXXFLAGS += -no_compact_unwind -std=c++11 -I$(MMG_HOME)/include -I$(PARMMG_HOME)/include -I$(PARMETIS_HOME)/include -I$(MPICH_HOME)/include -I$(HDF5_HOME)/include -I$(METIS_HOME)/include -I$(XML_HOME) -I$(BOOST_HOME) -I$(CGNS_HOME)

LDFLAGS += -L$(MMG_HOME)/lib -L$(PARMMG_HOME)/lib -L$(PARMETIS_HOME)/lib -L$(METIS_HOME)/lib -L$(MPICH_HOME)/lib -L$(HDF5_HOME)/lib -L$(XML_HOME) -L$(BOOST_HOME)/stage/lib
LDLIBS += -lmetis -lparmetis -lhdf5 -lmpi -llapack -lblas -lmmg -lparmmg

CC = /Users/dekelsch/Software/mpich-3.3.1/mpich-3.1.1-install/bin/mpic++

