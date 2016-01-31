#SYSTEM = LINUX
#SYSTEM = OSX
#SYSTEM = VSC3
#SYSTEM = MINGW

include ../Makefile.system

FFTW_BASE = $(HOME)
HDF5_BASE = $(HOME)
FC = gfortran 
FOPTS = -fopenmp -O3 -fPIC
CPP = cpp-4.9 -std=c89  -D_OPENMP
MAKE_SHARED_LIB = gfortran -fopenmp -shared 


TSSM_LIBS = -L. -ltssm
TSSM = libtssm.a
TSSM_SHARED_LIB = libtssm.so

ifeq ($(SYSTEM),VSC3)
    FFTW_BASE = /opt/sw/x86_64/glibc-2.12/ivybridge-ep/fftw/3.3.4/MPI_DP/intel-14.0.2
    HDF5_BASE = /opt/sw/x86_64/glibc-2.12/ivybridge-ep/hdf5/1.8.14/intel-14.0.2
    FC = ifort
    FOPTS =  -openmp -O3
    TSSM = libtssm.a
    TSSM_SHARED_LIB = libtssm.so
endif

ifeq ($(SYSTEM),OSX) 
    MAKE_SHARED_LIB = gfortran -fopenmp -dynamiclib
    CPP = cpp -std=c89  -D_OPENMP
    TSSM_SHARED_LIB = libtssm.dylib
all: shared_lib
endif   

ifeq ($(SYSTEM),MINGW)  
    FOPTS = -fopenmp -O3 
    TSSM_SHARED_LIB = libtssm.dll
    CPP = /c/MinGW/bin/cpp.exe -std=c89 -D_OPENMP
all: shared_lib
endif

ifeq ($(SYSTEM),LINUX) 
all: shared_lib
endif   


FFTW_INCLUDE = $(FFTW_BASE)/include 
FFTW_LIB = $(FFTW_BASE)/lib
FFTW_LIBS = -L$(FFTW_LIB) -lfftw3_omp -lfftw3 

HDF5_INCLUDE = $(HDF5_BASE)/include
HDF5_LIB = $(HDF5_BASE)/lib
HDF5_LIBS = -L$(HDF5_LIB) -lhdf5_fortran -lhdf5

ifeq ($(JULIA_BUILD),1)
    FOPTS = -fopenmp -O3 -fPIC -D_NO_HDF5_ 
    CPP = cpp -std=c89 -D_NO_HDF5_ -D_OPENMP 
    FFTW_INCLUDE = ../../include 
    HDF5_INCLUDE = ../../include 
    HDF5_LIBS =  
    FFTW_LIBS = /usr/lib/x86_64-linux-gnu/libfftw3_omp.so.3 /usr/lib/x86_64-linux-gnu/libfftw3.so.3 
else
    ifneq ($(SYSTEM),MINGW)  
        OBJS += tssm_hdf5.o 
    endif    
endif


include ../Makefile.common

shared_lib: $(TSSM_SHARED_LIB)

$(TSSM_SHARED_LIB):  $(OBJS)
	$(MAKE_SHARED_LIB) -o $(TSSM_SHARED_LIB) $(OBJS) $(FFTW_LIBS) $(HDF5_LIBS) 















