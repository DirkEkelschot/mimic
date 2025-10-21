.PHONY: default
default: build test

PYTHON = /Users/dekelsch/brew-4.1.24/Cellar/python@3.11/3.11.10/bin/python3.11
NUMPY_INCLUDE  = ${shell ${PYTHON} -c 'import numpy; print( numpy.get_include() )'}
MPI4PY_INCLUDE = ${shell ${PYTHON} -c 'import mpi4py; print( mpi4py.get_include() )'}
PYTHON_INCLUDE = /Users/dekelsch/brew-4.1.24/Cellar/python@3.11/3.11.10/Frameworks/Python.framework/Versions/3.11/include/python3.11
PYTHON_LIB = /Users/dekelsch/brew-4.1.24/Cellar/python@3.11/3.11.10/Frameworks/Python.framework/Versions/3.11/lib

PARMETIS_HOME = /Users/dekelsch/Software/parmetis-install
METIS_HOME = /Users/dekelsch/Software/metis-install
BOOST_HOME = /Users/dekelsch/Software/boost_1_83_0/boost-install
XML_HOME = /home1/dekelsch/Software/tinyxml
HDF5_HOME = /Users/dekelsch/Software/hdf5-1.14.2/install
ThirdParty_HOME = /Users/dekelsch/mimic-newmaster/build/ThirdParty/dist

CXXFLAGS += -fPIC -std=c++17 -g -DMPI_NO_CPPBIND -I$(PARMETIS_HOME)/include -I$(METIS_HOME)/include -I$(BOOST_HOME)/include -I$(XML_HOME) -I${HDF5_HOME}/include -I${ThirdParty_HOME}/include

LDFLAGS += -L$(PARMETIS_HOME)/lib -L$(METIS_HOME)/lib -L$(BOOST_HOME)/stage/lib -L$(XML_HOME) -L${HDF5_HOME}/lib -L${ThirdParty_HOME}/lib

$(info ${ThirdParty_HOME}/include)

LDLIBS += -lmpi -lparmetis -lmetis -lhdf5 -mkl -lmmg -lparmmg

CC = mpicc

MPICC = mpic++
BREW = $(HOME)/brew-4.1.24/Cellar
LDFLAGS += -L${PYTHON_LIB} -shared -ldl -framework CoreFoundation -lpython3.11m 
SOURCES=energy.cpp ../src/*.cpp
.PHONY: build 
build: mimic.so
mimic.so: mimic.cpp
	${MPICC} ${CXXFLAGS} -I${NUMPY_INCLUDE} -I${MPI4PY_INCLUDE} -I${PYTHON_INCLUDE} -o $@ $< $(SOURCES) ${LDFLAGS} ${LDLIBS}


MPIEXEC = mpiexec
NP_FLAG = -n
NP = 4
.PHONY: test
test: build
	${MPIEXEC} -f hostfile ${NP_FLAG} ${NP} ${PYTHON} mimic.py


.PHONY: clean
clean:
	${RM} *.png*
	${RM} *.dat*
	${RM} mimic.so
