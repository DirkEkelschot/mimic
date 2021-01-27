SRC 	 = src
BIN 	 = bin
TESGRAD  = tests/test_grad
OBJ      = $(SRC)/*.cpp
MAIN_OBJ = main.cpp
TES_OBJ  = $(TESGRAD)/main_test.cpp

include module.mk

all:	install tests
	
install:makebin
	$(CC) $(CXXFLAGS) $(OBJ) $(MAIN_OBJ) -o $(BIN)/adapt $(LDFLAGS) $(LDLIBS)
	cp -r $(BIN)/adapt .

makebin:
	mkdir -p $(BIN)
tests:
	make -C tests/test_grad
clean:	
	rm -rf adapt $(SRC)/*.o $(SRC)/*.mod $(BIN) $(TESGRAD)/testing
