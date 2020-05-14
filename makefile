OBJECTS = adapt_output.cpp adapt_compute.cpp adapt_partition.cpp adapt_operations.cpp main.cpp

PARMETIS_HOME = /Users/dekelschot/Software/parmetis-4.0.3/installation
METIS_HOME = /Users/dekelschot/Software/parmetis-4.0.3/metis/installation
HDF5_HOME = /Users/dekelschot/Software/hdf5-1.10.6/hdf5
MPICH_HOME = /usr/local/

CXXFLAGS += -std=c++11 -I$(PARMETIS_HOME)/include -I$(MPICH_HOME)/include -I$(HDF5_HOME)/include -I$(METIS_HOME)/include

LDFLAGS += -L$(PARMETIS_HOME)/lib -L$(METIS_HOME)/lib -L$(MPICH_HOME)/lib -L$(HDF5_HOME)/lib

LDLIBS += -lmetis -lparmetis -lhdf5 -lmpi -llapack -lblas

all:
	mpic++ $(CXXFLAGS) $(OBJECTS) -o adapt $(LDFLAGS) $(LDLIBS)
#	rm -rf *.o *.mod
