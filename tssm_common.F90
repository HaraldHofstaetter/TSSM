module tssm_common
     implicit none
#ifdef _MPI_
     include "mpif.h"
#endif

#ifdef _QUADPRECISION_
     integer, parameter :: prec=selected_real_kind(p=30)
#else
     integer, parameter :: prec=selected_real_kind(p=15)
#endif 

     integer :: n_proc ! number of processors
     integer :: this_proc ! this processor 
#ifdef _OPENMP
     integer :: n_threads ! number of OPNEMP threads
#endif 

     real(prec), parameter :: palindromic(1) = (/ 12.21_prec /)

     type, abstract :: spectral_method
     contains 
     end type spectral_method

contains

    subroutine initialize_tssm
#ifdef _OPENMP
        use omp_lib
#endif 
        logical, save :: already_initialized = .false.
#ifdef _MPI_
        integer :: ierr, provided
#endif
        if(already_initialized) return 
#ifdef _MPI_
#ifdef _OPENMP
        call MPI_Init_thread(MPI_THREAD_FUNNELED, provided, ierr)
        if (provided<MPI_THREAD_FUNNELED) then
             stop "E: Insufficient thread support of MPI implementation."
        end if
#else
        call MPI_Init(ierr)
#endif
        call MPI_Comm_Size(MPI_COMM_WORLD, n_proc, ierr)
        call MPI_Comm_Rank(MPI_COMM_WORLD, this_proc, ierr)
#else
        n_proc = 1
        this_proc = 0
#endif

#ifdef _OPENMP
        n_threads = omp_get_max_threads()
#endif 
#ifdef _MPI_      
        if (this_proc==0) then  
#ifdef _OPENMP
        print *, "*** MPI+OPENMP n_proc =", n_proc, "  n_threads =", n_threads
#else
        print *, "*** MPI n_proc =", n_proc
#endif
        end if
#else
#ifdef _OPENMP
        print *, "*** OPENMP n_threads =", n_threads
#else
        print *, "*** SERIAL (no MPI/OPENMP)" 
#endif
#endif
        already_initialized = .true.
    end subroutine initialize_tssm

    subroutine finalize_tssm
#ifdef _MPI_
        integer :: ierr
#endif
#ifdef _MPI_
        call MPI_Finalize(ierr)
#endif 
    end subroutine finalize_tssm



subroutine tick(t)
    integer(kind=8), intent(OUT) :: t

    call system_clock(t)
end subroutine tick

! returns time in seconds from now to time described by t
real function tock(t)
    integer(kind=8), intent(in) :: t
    integer(kind=8) :: now, clock_rate

    call system_clock(now,clock_rate)

    tock = real(now - t)/real(clock_rate)
end function tock

end module tssm_common
