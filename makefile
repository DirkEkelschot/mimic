OBJECTS = main.cpp

#HDF5_HOME = /u/smurman/share/eddy/hdf5-1.10.1
#PARMETIS_HOME = /nasa/modulefiles/sles12/parmetis/4.0.3
#MPI_HOME = /nasa/hpe/mpt/2.17r13
#CXXFLAGS += -std=c++11 -DMPI_NO_CPPBIND -I$(MPI_HOME)/include
CXXFLAGS += -std=c++11 -DMPI_NO_CPPBIND
#LDFLAGS += -L$(MPI_HOME)/lib

LDLIBS += -lmpi -lparmetis -lmetis -lhdf5 -mkl

all:
	icpc $(CXXFLAGS) main.cpp -o adapt $(LDLIBS)
#	rm -rf *.o *.mod
