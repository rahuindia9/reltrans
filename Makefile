BUILD  := build
ROOTDIR := .
HEADAS_LIB := ${HEADAS}/lib
HEADAS_INCLUDE := ${HEADAS}/include
DEBUG = 0
SANITIZE = 0

# The name to use for the reltrans library (must be different from reltrans, as
# libreltrans.so is the compiled reltrans library)
XSPEC_RELTRANS_NAME = xsreltrans

# This are configurable from the command line
TARGET = $(shell uname)
FC = gfortran

CFLAGS := -fno-omit-frame-pointer

FFLAGS := -DHAVE_INLINE -fPIC -fno-automatic -fno-second-underscore \
		  -fno-omit-frame-pointer \
		  -fopenmp \
		  -I$(BUILD)/include \
		  -I$(HEADAS_INCLUDE) \
		  -I$(HEADAS_INCLUDE)/fftw \
		  -J$(BUILD)/cache \
		  -I$(BUILD)/cache

LDFLAGS := -lXSFunctions -lXSModel -lfftw3 -lcfitsio \
		   -L$(BUILD)/lib -L$(HEADAS_LIB)

ifeq ($(DEBUG),1)
	FFLAGS += -g
	CFLAGS += -g
else
	FFLAGS += -O3
	CFLAGS += -O3
endif

ifeq ($(SANITIZE),1)
	FFLAGS += -fsanitize=address
	CFLAGS += -fsanitize=address
endif

ifeq ($(TARGET),Linux)
	FFLAGS += -shared -export-dynamic
	LDFLAGS += -lm -lpthread
	SHARED_EXT := so
	SED_INPLACE = sed -i
else
ifeq ($(TARGET),Darwin)
	FFLAGS += -dynamiclib
	LDFLAGS += -lgfortran
	SHARED_EXT := dylib
	# MacOS sed needs an extra useless argument
	SED_INPLACE = sed -i ''
endif
endif

# the path to the reltrans library for the -L linker flag
LIB_PATH := $(shell realpath $(BUILD))/lib
RELTRANS_SHARED_LIBRARY := $(BUILD)/lib/libreltrans.$(SHARED_EXT)

all: $(BUILD) $(RELTRANS_SHARED_LIBRARY)

exe: $(BUILD)/bin/relcli

$(BUILD)/bin/relcli: ./utils/cli.c $(BUILD)/lib/libreltrans.$(SHARED_EXT)
	$(CC) $(CFLAGS) utils/cli.c -o $(BUILD)/bin/relcli \
		-L$(BUILD)/lib -lgfortran -lc -lm -lmvec \
		-Wl,-rpath,'$$ORIGIN/../lib' -lreltrans

$(RELTRANS_SHARED_LIBRARY): $(BUILD) $(BUILD)/cache/wrappers.o
	$(FC) $(FFLAGS) $(BUILD)/cache/wrappers.o -o $@ $(LDFLAGS)

$(BUILD)/cache/%.o: $(ROOTDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)/bin $(BUILD)/lib $(BUILD)/include $(BUILD)/share $(BUILD)/cache

.PHONY: clean
clean:
	rm -rf $(BUILD)

.PHONY: format
format:
	clang-format -i ./utils/cli.c

.PHONY: xspec
xspec: $(RELTRANS_SHARED_LIBRARY) xspec/lmodel_reltrans.dat xspec/compile_reltrans.xcm
	# Copy the necessary XSPEC files into the build directory
	cp -r xspec $(BUILD)/xspec
	cd $(BUILD)/xspec && \
		echo "initpackage $(XSPEC_RELTRANS_NAME) lmodel_reltrans.dat .\n exit" | xspec
	# Delete the immediately compiled shared library, so we can recompile it
	# with our options.
	rm $(BUILD)/xspec/libxsreltrans.$(SHARED_EXT)
	# Patch the XSPEC generated Makefile so that it uses the shared library
	# compiled outside of XSPEC
	$(SED_INPLACE) 's|-lXSFunctions|-lXSFunctions -L$(LIB_PATH) -Wl,-rpath,"$(LIB_PATH)" -lreltrans|g' \
		$(BUILD)/xspec/Makefile
	# Set the library name
	$(SED_INPLACE) 's|{LIBRARY_NAME}|$(XSPEC_RELTRANS_NAME)|g' $(BUILD)/xspec/compile_reltrans.xcm
	# Compile and pray XSPEC is happy
	cd $(BUILD)/xspec && xspec - compile_reltrans.xcm

.PHONY: tables
tables:
	# Normalise the tables
	python3 ./renormalise_table.py
