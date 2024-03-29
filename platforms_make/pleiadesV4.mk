PARMMG_HOME = /home1/dekelsch/Software/parmmg_1.4/build
MMG_HOME = /home1/dekelsch/Software/parmmg_1.4/build/Mmg-prefix/src/Mmg-build
PARMETIS_HOME = /nobackupp19/dekelsch/parmetis-4.0.3/parmetis-install
METIS_HOME = /nobackupp19/dekelsch/parmetis-4.0.3/metis/metis-install
CGNS_HOME = /nobackupp19/dekelsch/CGNS-3.3.0/src/install
BOOST_HOME = /home1/dekelsch/Software/boost_1_71_0
XML_HOME = /home1/dekelsch/Software/tinyxml
CXXFLAGS += -std=c++11 -g -DMPI_NO_CPPBIND -I$(MMG_HOME)/include -I$(PARMMG_HOME)/include -I$(PARMETIS_HOME)/include -I$(METIS_HOME)/include -I$(BOOST_HOME) -I$(XML_HOME) -I$(CGNS_HOME)/include
LDFLAGS += -L$(MMG_HOME)/lib -L$(PARMMG_HOME)/lib -L$(PARMETIS_HOME)/lib -L$(METIS_HOME)/lib -L$(BOOST_HOME)/stage/lib -L$(XML_HOME) -L$(CGNS_HOME)/lib

LDLIBS += -lmpi -lparmetis -lmetis -lhdf5 -mkl -lparmmg -lmmg -lcgns

CC = icpc
