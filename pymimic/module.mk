PARMETIS_HOME = /nobackupp19/dekelsch/parmetis-4.0.3/parmetis-install
METIS_HOME = /nobackupp19/dekelsch/parmetis-4.0.3/metis/metis-install
BOOST_HOME = /nobackupp19/dekelsch/boost_1_71_0/boost-install
XML_HOME = /home1/dekelsch/Software/tinyxml
CXXFLAGS += -std=c++17 -g -DMPI_NO_CPPBIND -I$(PARMETIS_HOME)/include -I$(METIS_HOME)/include -I$(BOOST_HOME)/include -I$(XML_HOME)
LDFLAGS += -L$(PARMETIS_HOME)/lib -L$(METIS_HOME)/lib -L$(BOOST_HOME)/stage/lib -L$(XML_HOME}

LDLIBS += -lmpi -lparmetis -lmetis -lhdf5 -mkl -lparmmg -lmmg -lcgns

CC = icpc
