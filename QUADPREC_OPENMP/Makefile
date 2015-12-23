#SYSTEM = LINUX
#SYSTEM = OSX
#SYSTEM = VSC3
#SYSTEM = MINGW

include ../Makefile.system

FFTW_BASE = $(HOME)
HDF5_BASE = $(HOME)
FC = gfortran 
FOPTS = -fopenmp -O3 -fPIC -D_QUADPRECISION_
CPP = cpp-4.9 -std=c89  -D_OPENMP -D_QUADPRECISION_
MAKE_SHARED_LIB = gfortran -fopenmp -shared 


TSSM_LIBS = -L. -ltssmq
TSSM = libtssmq.a
TSSM_SHARED_LIB = libtssmq.so

ifeq ($(SYSTEM),VSC3)
    FFTW_BASE = /opt/sw/x86_64/glibc-2.12/ivybridge-ep/fftw/3.3.4/MPI_DP/intel-14.0.2
    HDF5_BASE = /opt/sw/x86_64/glibc-2.12/ivybridge-ep/hdf5/1.8.14/intel-14.0.2
    FC = ifort
    FOPTS =  -openmp -O3 -D_QUADPRECISION_
    TSSM = libtssmq.a
    TSSM_SHARED_LIB = libtssmq.so
endif

ifeq ($(SYSTEM),OSX) 
    MAKE_SHARED_LIB = gfortran -fopenmp -dynamiclib
    CPP = cpp -std=c89  -D_OPENMP -D_QUADPRECISION_
    TSSM_SHARED_LIB = libtssmq.dylib
all: shared_lib
endif   

ifeq ($(SYSTEM),MINGW)  
    FOPTS = -fopenmp -O3 -D_QUADPRECISION_
    TSSM_SHARED_LIB = libtssmq.dll
    CPP = /c/MinGW/bin/cpp.exe -std=c89 -D_OPENMP
all: shared_lib
endif

ifeq ($(SYSTEM),LINUX) 
all: shared_lib
endif   


FFTW_INCLUDE = $(FFTW_BASE)/include 
FFTW_LIB = $(FFTW_BASE)/lib
FFTW_LIBS = -L$(FFTW_LIB) -lfftw3q_omp -lfftw3_omp -lfftw3q -lfftw3 -lquadmath -lm

HDF5_INCLUDE = $(HDF5_BASE)/include
HDF5_LIB = $(HDF5_BASE)/lib
HDF5_LIBS = -L$(HDF5_LIB) -lhdf5_fortran -lhdf5

include ../Makefile.common

shared_lib: $(TSSM_SHARED_LIB)

$(TSSM_SHARED_LIB):  $(OBJS)
	$(MAKE_SHARED_LIB) -o $(TSSM_SHARED_LIB) $(OBJS) $(FFTW_LIBS) $(HDF5_LIBS) 















