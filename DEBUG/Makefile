SYSTEM = LINUX
#SYSTEM = OSX
#SYSTEM = VSC3
#SYSTEM = MINGW

#include ../Makefile.system

FFTW_BASE = $(HOME)
HDF5_BASE = $(HOME)
FC = gfortran 
FOPTS = -g -fPIC
CPP = cpp -std=c89  
MAKE_SHARED_LIB = gfortran -shared 


TSSM_LIBS = -L. -ltssm
TSSM = libtssm.a
TSSM_SHARED_LIB = libtssm.so

ifeq ($(SYSTEM),VSC3)
    FFTW_BASE = /opt/sw/x86_64/glibc-2.12/ivybridge-ep/fftw/3.3.4/MPI_DP/intel-14.0.2
    HDF5_BASE = /opt/sw/x86_64/glibc-2.12/ivybridge-ep/hdf5/1.8.14/intel-14.0.2
    FC = ifort
    FOPTS =  -g -fPIC
    TSSM = libtssm.a
    TSSM_SHARED_LIB = libtssm.so
    MAKE_SHARED_LIB = ifort -shared 
all: shared_lib
endif

ifeq ($(SYSTEM),OSX) 
    MAKE_SHARED_LIB = gfortran -dynamiclib
    CPP = cpp -std=c89
    TSSM_SHARED_LIB = libtssm.dylib
all: shared_lib
endif   

ifeq ($(SYSTEM),MINGW)  
    FOPTS = -g 
    TSSM_SHARED_LIB = libtssm.dll
    CPP = /c/MinGW/bin/cpp.exe -std=c89
all: shared_lib
endif

ifeq ($(SYSTEM),LINUX) 
    CPP = cpp -std=c89
all: shared_lib
endif   


FFTW_INCLUDE = $(FFTW_BASE)/include 
FFTW_LIB = $(FFTW_BASE)/lib
FFTW_LIBS = -L$(FFTW_LIB)  -lfftw3 

HDF5_INCLUDE = $(HDF5_BASE)/include
HDF5_LIB = $(HDF5_BASE)/lib
HDF5_LIBS = -L$(HDF5_LIB) -lhdf5_fortran -lhdf5

ifeq ($(SYSTEM),LINUX) 
    FFTW_INCLUDE = ../../include 
    HDF5_INCLUDE = ../../include 
    #HDF5_LIBS =  $(HDF5_LIB)/libhdf5_fortran.so \
    #             $(HDF5_LIB)/libhdf5.so
    #FFTW_LIBS = $(FFTW_LIB)/libfftw3.so.3 $(FFTW_LIB)/libfftw3_omp.so.3 
    HDF5_LIBS =  /usr/lib/x86_64-linux-gnu/libhdf5_serial_fortran.so.10 \
                 /usr/lib/x86_64-linux-gnu/libhdf5_serial.so.10
    FFTW_LIBS = /usr/lib/x86_64-linux-gnu/libfftw3_omp.so.3 /usr/lib/x86_64-linux-gnu/libfftw3.so.3 
endif


ifeq ($(JULIABOX_BUILD),1) # settings for the JuliaBox
    FOPTS = -g -fPIC 
    CPP = cpp -std=c89
    FFTW_INCLUDE = ../../include 
    HDF5_INCLUDE = ../../include 
    HDF5_LIBS =  /usr/lib/x86_64-linux-gnu/libhdf5_fortran.so.7 \
                 /usr/lib/x86_64-linux-gnu/libhdf5.so.7
    FFTW_LIBS =  /usr/lib/x86_64-linux-gnu/libfftw3.so.3 
endif



include ../Makefile.common

shared_lib: $(TSSM_SHARED_LIB)

$(TSSM_SHARED_LIB):  $(OBJS)
	$(MAKE_SHARED_LIB) -o $(TSSM_SHARED_LIB) $(OBJS) $(FFTW_LIBS) $(HDF5_LIBS) 















