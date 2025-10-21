
BUILD_DIR ?= build

BUILD_FLAGS=
INSTDM =  777
INSTFM =  777

ifndef US3D_HOME
$(error US3D_HOME is not set)
endif

ifndef PREFIX
ifneq ($(US3D_PLUGINS_DIR),)
PREFIX := $(US3D_PLUGINS_DIR)
else
PREFIX := $(US3D_HOME)/plugins
endif
endif

ifeq ($(DEBUG),true)
BUILD_FLAGS += --debug
endif

ifeq ($(DEBUG),on)
OPTIONS += --debug
#OPTIONS += --fcflags=-DUS3D_DEBUGGING
 OPTIONS += --fcflags="-DUS3D_DEBUGGING -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -Wunused-dummy-argument -Wno-unused-variable -Wno-array-temporaries -Wuninitialized -I$(ICARUS_HOME)/include/modfiles -ffree-form"
endif
#OPTIONS += --fcflags="-no-wrap-margin"

OPTIONS += --fcflags="-no-wrap-margin"

OBJECTS =   \
        src/compute_error_indicator_global.F90 \
	src/compute_error_indicator_io.F90 \
	src/compute_error_indicator_main.F90 \
	src/compute_error_indicator_driver.F90 \
	src/compute_error_indicator_setup.F90 \
	us3d_user.F90

TARGETS = $(BUILD_DIR)/general/compute-error-indicator-general.so

.NOTPARALLEL:

.PHONY: all
all: $(TARGETS)

$(BUILD_DIR)/general/compute-error-indicator-general.so: $(OBJECTS)
	#mkdir -p $(BUILD_DIR)/general
	#cd $(BUILD_DIR)/general
	$(US3D_HOME)/bin/us3d-build-user $(OPTIONS) \
	--output=compute-error-indicator-general.so $(OBJECTS);
	rm *.mod

install:
	install -m $(INSTDM) -d $(PREFIX)/compute-error-indicator
	install -m $(INSTDM) -d $(PREFIX)/compute-error-indicator/lib
	install -m $(INSTDM) -d $(PREFIX)/compute-error-indicator/inp
	install -m $(INSTFM) compute-error-indicator-general.so $(PREFIX)/compute-error-indicator/lib

clean:
	if [ -d "$(BUILD_DIR)" ] && [ -n "$(BUILD_DIR)" ] ; then \
           /bin/rm -rf "$(BUILD_DIR)"; \
        fi
	rm *.mod
	rm *.o

.PHONY: all install clean distclean
