
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
#OPTIONS += --fcflags="-no-wrap-margin"
OPTIONS += --fcflags=" -Werror"
# OPTIONS += --fcflags=" -I/$(SOFT_HOME)/metis/include/ -I/$(SOFT_HOME)/hdf5/include/"
# OPTIONS += --ldflags=" -L$(SOFT_HOME)/metis/lib/ -lparmetis -L/$(SOFT_HOME)/hdf5/lib/ -lhdf5"

OBJECTS =  us3d_user.F90

TARGETS = libus3d_user-general.so

.NOTPARALLEL:

.DEFAULT_GOAL := all
.PHONY: all
all: $(TARGETS)

# Plugins should define this target to create an info file that will help with
# the build process.  By making and parsing this file, we can determine some
# additional information about the plugin name, version, etc.  The plugin-info
# file should either be dynamically generated and not stored as a part of the
# repository itself, or it should be manually updated when things like the
# version number change.

# .PHONY: libus3d_user-general-general.so
# libus3d_user-general.so: $(OBJECTS)
# 	$(US3D_HOME)/bin/us3d-build-user \
# 	--output=libus3d_user-general-general.so $(OBJECTS);
# 	rm *.mod

.PHONY: libus3d_user-general-general.so
libus3d_user-general.so: $(OBJECTS)
	mpicxx -I$(SOFT_HOME)/metis/include -I$(SOFT_HOME)/hdf5/include -c ../adapt_compute.cpp -o adapt_compute.o
	$(US3D_HOME)/bin/us3d-build-user --ldflags="adapt_compute.o -lstdc++" \
	--output=libus3d_user-general.so $(OBJECTS);
	rm *.mod

.PHONY: all install clean distclean
