This test tests the METRIC reconstruction in parallel. For now the data is reduced to the vertices which for now is still dependent on the number of procs. In order to reproduce the reference data run:
mpiexec -np 4 ./adapt ../bin/test3
