OBJ = $(SRC)/*.cpp

SRC = src
BIN = bin
MAIN_OBJ = main.cpp
MAIN_API = main_api.cpp
EXEC = adapt
TESTBIN = tests/bin
TEST1 = tests/test1
TEST2 = tests/test2
TEST3 = tests/test3

include module.mk

all:		install
	
install:	makebin
		$(CC) $(CXXFLAGS) $(OBJ) $(MAIN_OBJ) -o $(BIN)/adapt $(LDFLAGS) $(LDLIBS)
		cp -r $(BIN)/$(EXEC) .

makebin:
		mkdir -p $(BIN)
test:
		make -C $(TEST1)
		make -C $(TEST2)
		make -C $(TEST3)
clean:
		rm -rf $(EXEC) *.dat $(SRC)/*.o $(SRC)/*.mod $(BIN) $(TESTBIN) grid_madam.h5
lib:
		rm -rf *.o *.a
		$(CC) $(CXXFLAGS) -fPIC -c $(SRC)/*.cpp main_api.cpp
		ar -rc libmadam.a \
		adapt_bltopology.o \
		adapt_boundary.o \
		adapt_compute.o \
		adapt_geometry.o \
		adapt_io.o \
		adapt_math.o \
		adapt_operations.o \
		adapt_output.o \
		adapt_parops.o \
		adapt_partition.o \
		adapt_recongrad.o \
		adapt_schedule.o \
		adapt_topology.o \
		hex2tet.o \
		main_api.o
