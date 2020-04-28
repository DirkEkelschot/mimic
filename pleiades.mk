OBJECTS = us3d_partitioning.cpp us3d_ops.cpp main.cpp

HDF5_HOME = /u/smurman/share/eddy/hdf5-1.10.1
PARMETIS_HOME = /u/smurman/share/eddy/parmetis-4.0.3
MPI_HOME = /nasa/hpe/mpt/2.17r13
CXXFLAGS += -std=c++11 -DMPI_NO_CPPBIND -I$(MPI_HOME)/include

LDFLAGS += -L$(MPI_HOME)/lib

LDLIBS += -lmpi -lparmetis -lmetis -lhdf5 -mkl

all:
	icpc $(CXXFLAGS) $(OBJECTS) -o adapt $(LDFLAGS) $(LDLIBS)
#	rm -rf *.o *.mod
