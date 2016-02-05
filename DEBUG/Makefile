#SYSTEM = LINUX
#SYSTEM = OSX
#SYSTEM = VSC3
#SYSTEM = MINGW

include ../Makefile.system

FFTW_BASE = $(HOME)
HDF5_BASE = $(HOME)
FC = gfortran 
FOPTS = -fcheck=all -g -fvisibility=hidden -fPIC
CPP = cpp-4.9 -std=c89 
#CPP = gfortran -E 
MAKE_SHARED_LIB = gfortran -shared


TSSM_LIBS = -L. -ltssm
TSSM = libtssm.a
TSSM_SHARED_LIB = libtssm.so

ifeq ($(SYSTEM),VSC3)
    FFTW_BASE = /opt/sw/x86_64/glibc-2.12/ivybridge-ep/fftw/3.3.4/MPI_DP/intel-14.0.2
    HDF5_BASE = /opt/sw/x86_64/glibc-2.12/ivybridge-ep/hdf5/1.8.14/intel-14.0.2
    FC = ifort
    FOPTS = -g
    TSSM = libtssm.a
    TSSM_SHARED_LIB = libtssm.so
endif

ifeq ($(SYSTEM),OSX) 
    MAKE_SHARED_LIB = gfortran -dynamiclib
    CPP = cpp -std=c89 
    TSSM_SHARED_LIB = libtssm.dylib
all: shared_lib
endif   

ifeq ($(SYSTEM),MINGW)  
    FOPTS = -fcheck=all -g
    TSSM_SHARED_LIB = libtssm.dll
    CPP = /c/MinGW/bin/cpp.exe -std=c89
all: shared_lib
endif

ifeq ($(SYSTEM),LINUX) 
all: shared_lib
endif   


FFTW_INCLUDE = $(FFTW_BASE)/include 
FFTW_LIB = $(FFTW_BASE)/lib
#FFTW_LIBS = -L$(FFTW_LIB) -lfftw3 
FFTW_LIBS = $(FFTW_LIB)/libfftw3.so 

HDF5_INCLUDE = $(HDF5_BASE)/include
HDF5_LIB = $(HDF5_BASE)/lib
HDF5_LIBS = -L$(HDF5_LIB) -lhdf5_fortran -lhdf5

include ../Makefile.common

shared_lib: $(TSSM_SHARED_LIB)

$(TSSM_SHARED_LIB):  $(OBJS)
	$(MAKE_SHARED_LIB) -o $(TSSM_SHARED_LIB) $(OBJS) $(FFTW_LIBS) $(HDF5_LIBS) 



















