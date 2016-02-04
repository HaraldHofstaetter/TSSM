#SYSTEM = VSC3
include ../Makefile.system

FFTW_BASE = $(HOME)
HDF5_BASE = $(HOME)
FC = gfortran 
FOPTS = -fcheck=all -g -fPIC -D_QUADPRECISION_
CC = gcc
COPTS = -Wall -g -fPIC -D_QUADPRECISION_
CPP = cpp-4.9 -std=c89 -D_QUADPRECISION_

#CPP = gfortran -E -D_QUADPRECISION
MAKE_SHARED_LIB = gfortran -shared


TSSM_LIBS = -L. -ltssm
TSSM = libtssmq.a
TSSM_SHARED_LIB = libtssmq.so

ifeq ($(SYSTEM),VSC3)
    FFTW_BASE = /opt/sw/x86_64/glibc-2.12/ivybridge-ep/fftw/3.3.4/MPI_DP/intel-14.0.2
    HDF5_BASE = /opt/sw/x86_64/glibc-2.12/ivybridge-ep/hdf5/1.8.14/intel-14.0.2
    FC = ifort
    FOPTS = -g  -D_QUADPRECISION_
    CC = icc
    COPTS = -g  -D_QUADPRECISION_
    TSSM = libtssmq.a
    TSSM_SHARED_LIB = libtssmq.so
endif

ifeq ($(SYSTEM),OSX) 
    MAKE_SHARED_LIB = gfortran -dynamiclib
    CPP = cpp -std=c89 -D_QUADPRECISION_
    TSSM_SHARED_LIB = libtssmq.dylib
all: shared_lib
endif   

ifeq ($(SYSTEM),MINGW)  
    FOPTS = -fcheck=all -g -D_QUADPRECISION_
    TSSM_SHARED_LIB = libtssmq.dll
    CPP = /c/MinGW/bin/cpp.exe -std=c89 -D_QUADPRECISION_
all: shared_lib
endif

ifeq ($(SYSTEM),LINUX) 
all: shared_lib
endif   


FFTW_INCLUDE = $(FFTW_BASE)/include 
FFTW_LIB = $(FFTW_BASE)/lib
FFTW_LIBS = -L$(FFTW_LIB) -lfftw3q -lfftw3 -lquadmath -lm

HDF5_INCLUDE = $(HDF5_BASE)/include
HDF5_LIB = $(HDF5_BASE)/lib
HDF5_LIBS = -L$(HDF5_LIB) -lhdf5_fortran -lhdf5

include ../Makefile.common

shared_lib: $(TSSM_SHARED_LIB)

$(TSSM_SHARED_LIB):  $(OBJS) $(OBJS_F128WRAPPER)
	$(MAKE_SHARED_LIB) -o $(TSSM_SHARED_LIB) $(OBJS) $(OBJS_F128WRAPPER) $(FFTW_LIBS) $(HDF5_LIBS) 























