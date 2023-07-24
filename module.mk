PARMETIS_HOME = /Users/dekelschot/Software/ParMETIS/parmetis-install
METIS_HOME = /Users/dekelschot/Software/METIS/metis-install
HDF5_HOME = /Users/dekelschot/Software/hdf5-1.12.0/hdf5-install
MPICH_HOME = /Users/dekelschot/Software/mpich-3.3.1/mpich-3.3.1-install
PARMMG_HOME = /Users/dekelschot/Software/ParMmg/build
MMG_HOME = /Users/dekelschot/Software/ParMmg/build/Mmg-prefix/src/Mmg-build
XML_HOME = /Users/dekelschot/Software/tinyxml
BOOST_HOME = /Users/dekelschot/Software/boost_1_71_0
CGNS_HOME = /Users/dekelschot/Software/CGNS-3.3.0/build
GKLIB_HOME = /Users/dekelschot/Software/GKlib/GKlib-install

CXXFLAGS += -std=c++17 -I$(MMG_HOME)/include -I$(PARMMG_HOME)/include -I$(PARMETIS_HOME)/include -I$(MPICH_HOME)/include -I$(HDF5_HOME)/include -I$(METIS_HOME)/include -I$(XML_HOME) -I$(BOOST_HOME) -I$(CGNS_HOME)/include -I$(GKLIB_HOME)/include

LDFLAGS += -L$(MMG_HOME)/lib -L$(PARMMG_HOME)/lib -L$(PARMETIS_HOME)/lib -L$(METIS_HOME)/lib -L$(MPICH_HOME)/lib -L$(HDF5_HOME)/lib -L$(XML_HOME) -L$(CGNS_HOME)/lib -L$(GKLIB_HOME)/lib

LDLIBS += -lmetis -lparmetis -lhdf5 -lmpi -llapack -lblas -lmmg -lparmmg -lcgns -lgklib



CC = /Users/dekelschot/Software/mpich-3.3.1/mpich-3.3.1-install/bin/mpic++

