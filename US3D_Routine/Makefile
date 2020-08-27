
BUILD_DIR ?= build

BUILD_FLAGS=
INSTDM =  777
INSTFM =  777

ifeq ($(DEBUG),true)
BUILD_FLAGS += --debug
endif

ifeq ($(DEBUG),on)
OPTIONS += --debug
OPTIONS += --fcflags="-g3"
endif
OPTIONS += --fcflags=" -Werror"

OBJECTS =  us3d_user.F90

TARGETS = libus3d_user-general.so

.NOTPARALLEL:

.DEFAULT_GOAL := all
.PHONY: all
all: $(TARGETS)

.PHONY: libus3d_user-general-general.so
libus3d_user-general.so: $(OBJECTS)
	mpicxx -I$(SOFT_HOME)/metis/include -I$(SOFT_HOME)/hdf5/include -c ../adapt_compute.cpp -o adapt_compute.o
	$(US3D_HOME)/bin/us3d-build-user --ldflags="adapt_compute.o" \
	--output=libus3d_user-general.so $(OBJECTS);
	rm *.mod

.PHONY: all install clean distclean
