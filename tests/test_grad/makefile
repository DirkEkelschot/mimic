SRC_OBJ = ../../src/*.cpp
TES_OBJ = main_test.cpp

include ../../module.mk

test:
	$(CC) $(CXXFLAGS) $(SRC_OBJ) $(TES_OBJ) -o testing $(LDFLAGS) $(LDLIBS)

clean:	
	rm -rf testing
