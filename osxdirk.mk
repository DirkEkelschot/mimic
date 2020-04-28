OBJECTS = us3d_partitioning.cpp us3d_ops.cpp main.cpp

PARMETIS_HOME = /Users/dekelsch/Software/parmetis-4.0.3/parmetis-install
METIS_HOME = /Users/dekelsch/Software/parmetis-4.0.3/metis/metis-install
HDF5_HOME = /Users/dekelsch/Software/hdf5-1.12.0/hdf5-install
MPICH_HOME = /Users/dekelsch/Software/mpich-3.3.2/mpich-install

CXXFLAGS += -std=c++11 -I$(PARMETIS_HOME)/include -I$(MPICH_HOME)/include -I$(HDF5_HOME)/include -I$(METIS_HOME)/include

LDFLAGS += -L$(PARMETIS_HOME)/lib -L$(METIS_HOME)/lib -L$(MPICH_HOME)/lib -L$(HDF5_HOME)/lib

LDLIBS += -lmetis -lparmetis -lhdf5 -lmpi -llapack -lblas

all:
	/Users/dekelsch/Software/mpich-3.3.2/mpich-install/bin/mpic++ $(CXXFLAGS) $(OBJECTS) -o adapt $(LDFLAGS) $(LDLIBS)
#	rm -rf *.o *.mod
