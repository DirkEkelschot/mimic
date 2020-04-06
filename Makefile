OBJECTS = main.cpp

CXXFLAGS += -std=c++11 -I/u/smurman/share/eddy/hdf5-1.10.1/include -I/u/smurman/share/eddy/parmetis-4.0.3/include -fPIC -lm -lpthread -ldl -lutil -O2 -Wall -axAVX,CORE-AVX2,CORE-AVX512 -qopt-report-phase=vec -qopt-report=5

LDFLAGS += -L$(HOME)/lib -L/u/smurman/share/eddy/hdf5-1.10.1/lib -L/u/smurman/share/eddy/parmetis-4.0.3/lib -mkl

LDLIBS += -lmpi -lparmetis -lmetis -lhdf5 -lhdf5_hl -ldl

all:
	h5pcc $(CXXFLAGS) main.cpp -o adapt $(LDFLAGS) $(LDLIBS)
#	rm -rf *.o *.mod
