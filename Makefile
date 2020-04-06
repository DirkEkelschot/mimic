OBJECTS = main.cpp

CXXFLAGS += -std=c++11 -I/Users/dekelsch/Software/parmetis-4.0.3/parmetis-install/include -I/Users/dekelsch/Software/mpich-3.3.2/mpich-install/include -I/Users/dekelsch/Software/hdf5-1.12.0/hdf5-install/include -I/Users/dekelsch/Software/parmetis-4.0.3/metis/metis-install/include

LDFLAGS += -L/Users/dekelsch/Software/parmetis-4.0.3/parmetis-install/lib -L/Users/dekelsch/Software/parmetis-4.0.3/metis/metis-install/lib -L/Users/dekelsch/Software/mpich-3.3.2/mpich-install/lib -L/Users/dekelsch/Software/hdf5-1.12.0/hdf5-install/lib -lmetis -lparmetis -lhdf5 -lmpi -llapack -lblas

all:
	/Users/dekelsch/Software/mpich-3.3.2/mpich-install/bin/mpic++ $(CXXFLAGS) main.cpp -o compute_jacobians $(LDFLAGS)
#	rm -rf *.o *.mod
