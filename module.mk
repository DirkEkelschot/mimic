PARMMG_HOME = /home1/dekelsch/Software/parmmg_0920/build
MMG_HOME = /home1/dekelsch/Software/mmg_0920/build
PARMETIS_HOME = /nobackupp2/dekelsch/parmetis-4.0.3/parmetis-install
METIS_HOME = /nobackupp2/dekelsch/parmetis-4.0.3/metis/metis-install
CXXFLAGS += -std=c++11 -DMPI_NO_CPPBIND -I$(MMG_HOME)/include -I$(PARMMG_HOME)/include -I$(PARMETIS_HOME)/include -I$(METIS_HOME)/include -O3 -axCORE-AVX512,CORE-AVX2 -xAVX 
LDFLAGS += -L$(MMG_HOME)/lib -L$(PARMMG_HOME)/lib -L$(PARMETIS_HOME)/lib -L$(METIS_HOME)/lib
LDLIBS += -lmpi -lparmetis -lmetis -lhdf5 -mkl -lparmmg -lmmg

CC = icpc
