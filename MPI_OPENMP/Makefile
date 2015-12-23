#SYSTEM = VSC3
include ../Makefile.system

ifeq ($(SYSTEM),VSC3)
    FFTW_BASE = /opt/sw/x86_64/glibc-2.12/ivybridge-ep/fftw/3.3.4/MPI_DP/intel-14.0.2
    HDF5_BASE = /opt/sw/x86_64/glibc-2.12/ivybridge-ep/hdf5/1.8.14/intel-14.0.2
    FC = mpiifort
    FOPTS =  -D_MPI_ -openmp -O3 
else
    FFTW_BASE = $(HOME)
    HDF5_BASE = $(HOME)
    FC = mpifort 
    FOPTS =  -D_MPI_ -fopenmp -O3 
endif

CPP = cpp -std=c89 -D_MPI_ -D_OPENMP

TSSM = libtssm_mpi_omp.a
TSSM_LIBS = -L. -ltssm_mpi_omp


FFTW_INCLUDE = $(FFTW_BASE)/include 
FFTW_LIB = $(FFTW_BASE)/lib
FFTW_LIBS = -L$(FFTW_LIB) -lfftw3_mpi -lfftw3_omp -lfftw3 

HDF5_INCLUDE = $(HDF5_BASE)/include
HDF5_LIB = $(HDF5_BASE)/lib
HDF5_LIBS = -L$(HDF5_LIB) -lhdf5_fortran -lhdf5

include ../Makefile.common
