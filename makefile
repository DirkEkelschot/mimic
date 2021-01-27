SRC 	 = src
BIN 	 = bin
TESGRAD  = tests/test_grad
OBJ      = $(SRC)/*.cpp
MAIN_OBJ = main.cpp
TES_OBJ  = $(TESGRAD)/main_test.cpp

include module.mk

all:	makebin build_app tests
	
makebin:
	mkdir -p $(BIN)

build_app:
	$(CC) $(CXXFLAGS) $(OBJ) $(MAIN_OBJ) -o $(BIN)/adapt $(LDFLAGS) $(LDLIBS)

tests:
	make -C tests/test_grad
clean:	
	rm -rf $(SRC)/*.o $(SRC)/*.mod $(BIN) $(TESGRAD)/testing
