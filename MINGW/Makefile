FFTW_BASE = $(HOME)/MinGW
HDF5_BASE = $(HOME)/MinGW
FC = i686-w64-mingw32-gfortran -fopenmp
FOPTS = -O3  -D_NO_HDF5_
CPP = cpp-4.9 -std=c89 -D_NO_HDF5_ -D_OPENMP
MAKE_SHARED_LIB =  i686-w64-mingw32-gfortran -shared -fopenmp

TSSM_LIBS = -L. -ltssm
TSSM = libtssm.a
TSSM_SHARED_LIB = libtssm.dll

all: shared_lib


FFTW_INCLUDE = $(FFTW_BASE)/include 
FFTW_LIB = $(FFTW_BASE)/lib
FFTW_LIBS = -L$(FFTW_LIB) -lfftw3

HDF5_INCLUDE = $(HDF5_BASE)/include
HDF5_LIB = $(HDF5_BASE)/lib
#HDF5_LIBS = -L$(HDF5_LIB) -lhdf5_fortran -lhdf5
HDF5_LIBS = 

include ../Makefile.common

shared_lib: $(TSSM_SHARED_LIB)

$(TSSM_SHARED_LIB):  $(OBJS)
	$(MAKE_SHARED_LIB) -o $(TSSM_SHARED_LIB) $(OBJS) $(FFTW_LIBS) $(HDF5_LIBS) 








