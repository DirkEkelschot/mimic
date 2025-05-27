.PHONY: default
default: build test clean

PYTHON = /home1/dekelsch/local/bin/python3.12
#PYTHON_CONFIG = ${PYTHON} python-config
MPI4PY_INCLUDE = ${shell ${PYTHON} -c 'import mpi4py; print( mpi4py.get_include() )'}
NUMPY_INCLUDE  = ${shell ${PYTHON} -c 'import numpy; print( numpy.get_include() )'}

PYTHON_INCLUDE = /home1/dekelsch/local/include/python3.12/ 
PYTHON_LIB = /home1/dekelsch/local/lib/python3.12/

include module.mk

MPICC    = icc
CXXFLAGS += -fPIC -lm -lpthread -ldl -lutil -L${PYTHON_LIB} -O2 -Wall -axAVX,CORE-AVX2,CORE-AVX512 -qopt-report-phase=vec -qopt-report=5 
LDFLAGS += -shared -lm -lpthread -ldl -lutil -L${PYTHON_LIB}
SOURCES=energy.cpp src/*.cpp
LDLIBS  += -lmpi
#LDFLAGS += ${shell ${PYTHON_CONFIG} --libs}
#SO = ${shell ${PYTHON_CONFIG} --extension-suffix}
.PHONY: build
build: PartiClass.so
PartiClass.so: PyPartiClass.cpp
	${MPICC} ${CXXFLAGS} -I${NUMPY_INCLUDE} -I${MPI4PY_INCLUDE} -I${PYTHON_INCLUDE} -o $@ $< $(SOURCES) ${LDFLAGS} ${LDLIBS}


MPIEXEC = mpiexec
NP_FLAG = -n
NP = 4
.PHONY: test
test: build
	${MPIEXEC} ${NP_FLAG} ${NP} ${PYTHON} run_partiClass.py


.PHONY: clean
clean:
	${RM} PartiClass.so
