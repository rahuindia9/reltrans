BUILD  := build
ROOTDIR := ..
HEADAS_LIB := ${HEADAS}/lib
FFTWDIR := $(ROOTDIR)/fftw-3.3.9

FC := gfortran

CFLAGS := -g -fno-omit-frame-pointer

FFLAGS := -DHAVE_INLINE -g -fPIC -fno-automatic -fno-second-underscore \
		  -fopenmp \
		  -O3 \
		  -I$(BUILD)/include \
		  -I$(FFTWDIR)/api \
		  -J$(BUILD)/cache \
		  -I$(BUILD)/cache

LDFLAGS := -lXSFunctions -lXSModel -lfftw3 -lm -lcfitsio -lpthread \
		   -L$(BUILD)/lib -L$(HEADAS_LIB) -shared -export-dynamic

all: $(BUILD) $(BUILD)/lib/libreltrans.so

exe: $(BUILD)/bin/relcli

$(BUILD)/bin/relcli: ./cli.c $(BUILD)/lib/libreltrans.so
	$(CC) $(CFLAGS) cli.c -o $(BUILD)/bin/relcli -Wl,-rpath,$(BUILD)/lib \
		-L$(BUILD)/lib -lgfortran -lc -lreltrans -lm -lmvec

$(BUILD)/lib/libreltrans.so: $(FFTWDIR)/api/fftw3.f03 $(BUILD)/cache/wrappers.o
	$(LD) $(LDFLAGS) $(BUILD)/cache/wrappers.o -o $@

$(BUILD)/cache/%.o: $(ROOTDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(FFTWDIR)/api/fftw3.f03: $(FFTWDIR)/api/fftw3.f03.in

$(FFTWDIR)/api/fftw3.f03.in:
	tar xf ./$(FFTWDIR).tar.gz
	( cd $(FFTWDIR) && ./configure && make -j $(shell nproc) )


$(BUILD):
	mkdir -p $(BUILD)/bin $(BUILD)/lib $(BUILD)/include $(BUILD)/share $(BUILD)/cache

.PHONY: clean
clean:
	rm -rf $(BUILD)

.PHONY: format
format:
	clang-format -i ./cli.c
