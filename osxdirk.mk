OBJECTS = adapt_output.cpp\
	  adapt_compute.cpp\
	  adapt_schedule.cpp\
	  adapt_operations.cpp\
	  hex2tet.cpp \
	  adapt_geometry.cpp \
	  adapt_math.cpp \
	  adapt_recongrad.cpp \
	  adapt_io.cpp \
 	  adapt_topology.cpp \
	  adapt_partition.cpp \
	  main.cpp

PARMETIS_HOME = /Users/dekelsch/Software/parmetis-4.0.3/parmetis-install
METIS_HOME = /Users/dekelsch/Software/parmetis-4.0.3/metis/metis-install
HDF5_HOME = /Users/dekelsch/Software/hdf5-1.12.0/hdf5-install
MPICH_HOME = /Users/dekelsch/Software/mpich-3.3.1/mpich-install
PARMMG_HOME = /Users/dekelsch/Software/parmmg/build
MMG_HOME = /Users/dekelsch/Software/parmmg/build/Mmg-prefix/src/Mmg-build
CXXFLAGS += -std=c++11 -I$(MMG_HOME)/include -I$(PARMMG_HOME)/include -I$(PARMETIS_HOME)/include -I$(MPICH_HOME)/include -I$(HDF5_HOME)/include -I$(METIS_HOME)/include

LDFLAGS += -L$(MMG_HOME)/lib -L$(PARMMG_HOME)/lib -L$(PARMETIS_HOME)/lib -L$(METIS_HOME)/lib -L$(MPICH_HOME)/lib -L$(HDF5_HOME)/lib

LDLIBS += -lmetis -lparmetis -lhdf5 -lmpi -llapack -lblas -lmmg -lparmmg

all:
	/Users/dekelsch/Software/mpich-3.3.1/mpich-install/bin/mpic++ $(CXXFLAGS) $(OBJECTS) -o adapt $(LDFLAGS) $(LDLIBS)
#	rm -rf *.o *.mod
