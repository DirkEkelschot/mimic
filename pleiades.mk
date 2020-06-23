OBJECTS = adapt_output.cpp\
	  adapt_compute.cpp\
	  adapt_part_func.cpp\
	  schedule.cpp\
	  adapt_operations.cpp\
	  main.cpp

CXXFLAGS += -std=c++11 -DMPI_NO_CPPBIND

LDLIBS += -lmpi -lparmetis -lmetis -lhdf5 -mkl

all:
	icpc $(CXXFLAGS) $(OBJECTS) -o adapt $(LDLIBS)
#	rm -rf *.o *.mod
