BUILD  := build
ROOTDIR := ..
HEADAS_LIB := ${HEADAS}/lib
HEADAS_INCLUDE := ${HEADAS}/include
DEBUG = 0
SANITIZE = 0

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
else
ifeq ($(TARGET),Darwin)
	FFLAGS += -dynamiclib
	LDFLAGS += -lgfortran
	SHARED_EXT := dylib
endif
endif

all: $(BUILD) $(BUILD)/lib/libreltrans.$(SHARED_EXT)

exe: $(BUILD)/bin/relcli

$(BUILD)/bin/relcli: ./cli.c $(BUILD)/lib/libreltrans.$(SHARED_EXT)
	$(CC) $(CFLAGS) cli.c -o $(BUILD)/bin/relcli -Wl,-rpath,$(BUILD)/lib \
		-L$(BUILD)/lib -lgfortran -lc -lreltrans -lm -lmvec

$(BUILD)/lib/libreltrans.$(SHARED_EXT): $(BUILD)/cache/wrappers.o
	$(FC) $(LDFLAGS) $(FFLAGS) $(BUILD)/cache/wrappers.o -o $@

$(BUILD)/cache/%.o: $(ROOTDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)/bin $(BUILD)/lib $(BUILD)/include $(BUILD)/share $(BUILD)/cache

.PHONY: clean
clean:
	rm -rf $(BUILD)

.PHONY: format
format:
	clang-format -i ./cli.c
