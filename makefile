OBJ = $(SRC)/*.cpp

SRC = src
BIN = bin
MAIN_OBJ = main.cpp
EXEC = adapt
TESTBIN = tests/bin
TEST1 = tests/test1
TEST2 = tests/test2

include module.mk

all:	install test
	
install:makebin
	$(CC) $(CXXFLAGS) $(OBJ) $(MAIN_OBJ) -o $(BIN)/adapt $(LDFLAGS) $(LDLIBS)
	cp -r $(BIN)/$(EXEC) .

makebin:
	mkdir -p $(BIN)
test:
	make -C $(TEST1)
	make -C $(TEST2)
clean:	
	rm -rf $(EXEC) *.dat $(SRC)/*.o $(SRC)/*.mod $(BIN) $(TESTBIN)
