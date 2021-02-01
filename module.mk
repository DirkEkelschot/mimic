PARMMG_HOME = /home1/dekelsch/Software/parmmg/build
MMG_HOME = /home1/dekelsch/Software/parmmg/build/Mmg-prefix/src/Mmg-build

CXXFLAGS += -std=c++11 -DMPI_NO_CPPBIND -I$(MMG_HOME)/include -I$(PARMMG_HOME)/include
LDFLAGS += -L$(MMG_HOME)/lib -L$(PARMMG_HOME)/lib
LDLIBS += -lmpi -lparmetis -lmetis -lhdf5 -mkl -lparmmg -lmmg

CC = icpc
