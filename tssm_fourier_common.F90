#ifdef _QUADPRECISION_
#define fftw_cleanup fftwq_cleanup
#define fftw_import_wisdom_from_filename fftwq_import_wisdom_from_filename
#define fftw_export_wisdom_to_filename  fftwq_export_wisdom_to_filename
#define fftw_init_threads fftwq_init_threads
#define fftw_plan_with_nthreads fftwq_plan_with_nthreads
#define fftw_cleanup_threads fftwq_cleanup_threads
#define fftw_destroy_plan fftwq_destroy_plan
#define fftw_free fftwq_free
#endif

#ifdef _QUADPRECISION_
module tssmq_fourier_common
    use tssmq_base
#else
module tssm_fourier_common
    use tssm_base
#endif    
    use, intrinsic :: iso_c_binding !needed for fftw3.f03
    implicit none
#if(defined(_MPI_))
    include 'fftw3-mpi.f03'
#elif (defined(_QUADPRECISION_))
    include 'fftw3.f03'
    include 'fftw3q.f03'
#else
    include 'fftw3.f03'
#endif 

    integer(C_INT) :: fftw_planning_rigor = FFTW_PATIENT
#ifdef _MPI_
    character(len=64, kind=c_char) :: fftw_wisdom_file = c_char_'fftw_mpi_wisdom' // c_null_char
#else 
#ifdef _QUADPRECISION_
        character(len=64, kind=c_char) :: fftw_wisdom_file = c_char_'fftwq_wisdom' // c_null_char
#else
        character(len=64, kind=c_char) :: fftw_wisdom_file = c_char_'fftw_wisdom' // c_null_char
#endif
#endif 
    
    ! boundary condition kinds 
    integer, parameter :: periodic = 0 ! standard discrete fourier transform (periodic boundary conditions)
    integer, parameter :: dirichlet = 1 ! discrete sine transform (homogeneous Dirichlet boundary conditions)
    integer, parameter :: neumann = 2 ! discrete cosine transform (homogeneous Neumann boundary conditions)

    real(kind=prec), parameter :: pi = 4.0_prec*atan(1.0_prec)

contains

    subroutine initialize_tssm_fourier
#ifdef _OPENMP    
        use omp_lib
        integer :: ierr
#endif
        logical, save :: already_initialized = .false.
        if(already_initialized) return 
        call initialize_tssm
#ifdef _OPENMP    
        ierr = fftw_init_threads()
#endif
#ifdef _MPI_
        call fftw_mpi_init
#endif
        call load_fftw_wisdom
#ifdef _OPENMP    
        call fftw_plan_with_nthreads(n_threads)
#endif
        already_initialized = .true.
    end subroutine initialize_tssm_fourier

    subroutine finalize_tssm_fourier
            call save_fftw_wisdom

#ifdef _OPENMP    
            call fftw_cleanup_threads
#endif            
#ifdef _MPI_
            call fftw_mpi_cleanup
#else 
            call fftw_cleanup
#endif    
            call finalize_tssm
    end subroutine finalize_tssm_fourier


    subroutine load_fftw_wisdom
        logical, save :: already_loaded = .false.
        integer(C_INT) :: ret

        if (already_loaded) return
#ifdef _MPI_
        if (this_proc==0) then
#endif
            ret = fftw_import_wisdom_from_filename(fftw_wisdom_file)
            if (ret .eq. 0) then
                print *, 'W: could not import fftw wisdom file'
            end if    
#ifdef _MPI_
        end if
        call fftw_mpi_broadcast_wisdom(MPI_COMM_WORLD)
#endif         
        already_loaded = .true.
    end subroutine load_fftw_wisdom

    subroutine save_fftw_wisdom
        integer(C_INT) :: ret
#ifdef _MPI_
        call fftw_mpi_gather_wisdom(MPI_COMM_WORLD)
        if (this_proc==0) then
#endif         
            ret = fftw_export_wisdom_to_filename(fftw_wisdom_file)
            if (ret .eq. 0) then
                print *, 'W: could not  export fftw wisdom file'
            end if    
#ifdef _MPI_
       end if
#endif         
    end subroutine save_fftw_wisdom

    subroutine get_eigenvalues(lambda, d, n, nmin, nmax, boundary_conditions)
        real(kind=prec), intent(out) :: lambda(nmin:nmax)
        real(kind=prec), intent(in) :: d
        integer, intent(in) :: n, nmin, nmax
        integer, optional, intent(in):: boundary_conditions 
        integer :: m,k

        if (present(boundary_conditions).and.((boundary_conditions==dirichlet) &
            .or.(boundary_conditions==neumann))) then
            lambda(nmin:nmax) = (/ (real(k, kind=prec)**2, k=nmin, nmax) /)
            lambda = -(pi/d)**2*lambda
        else
            m = min(n/2, nmax)
            if (m>=nmin) then
                lambda(nmin:m) = (/ (real(k, kind=prec)**2, k=nmin-1, m-1) /)
            end if    
            m = max(n/2+1, nmin)
            if (m<=nmax) then
                lambda(m:nmax) = (/ (real(k, kind=prec)**2, k=m-n-1, nmax-n-1) /)
            end if    
            lambda = (-(2.0_prec*pi/d)**2) * lambda
        endif
    end subroutine get_eigenvalues

    subroutine get_eigenvalues_d(lambda, d, n, nmin, nmax, boundary_conditions)
        real(kind=prec), intent(out) :: lambda(nmin:nmax)
        real(kind=prec), intent(in) :: d
        integer, intent(in) :: n, nmin, nmax
        integer, optional, intent(in):: boundary_conditions 
        integer :: m,k

        if (present(boundary_conditions).and.((boundary_conditions==dirichlet) &
            .or.(boundary_conditions==neumann))) then
            lambda(nmin:nmax) = (/ (real(k, kind=prec), k=nmin, nmax) /)
            lambda = (-pi/d)*lambda
        else
            m = min(n/2, nmax)
            if (m>=nmin) then
                lambda(nmin:m) = (/ (real(k, kind=prec), k=nmin-1, m-1) /)
            end if    
            m = max(n/2+1, nmin)
            if (m<=nmax) then
                lambda(m:nmax) = (/ (real(k, kind=prec), k=m-n-1, nmax-n-1) /)
            end if    
            lambda = (-2.0_prec*pi/d) * lambda
        endif
    end subroutine get_eigenvalues_d




#ifdef _MPI_
    subroutine get_scrambled_eigenvalues(lambda, d, n, nmin, nmax, alloc_size)
        real(kind=prec), intent(out) :: lambda(nmin:nmax)
        real(kind=prec), intent(in) :: d
        integer, intent(in) :: n, nmin, nmax
        integer(kind=C_SIZE_T), intent(in) :: alloc_size
        integer :: m,k

        type(C_PTR) :: p
        complex(kind=8), pointer  :: u(:)
        type(c_ptr) :: plan_forward_scrambled 
        type(c_ptr) :: plan_backward

        p = fftw_alloc_complex(alloc_size)
        call c_f_pointer(p, u, [alloc_size])

        plan_backward = fftw_mpi_plan_dft_1d ( int(n, kind=C_SIZE_T), u, u, &
               MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE )

        plan_forward_scrambled = fftw_mpi_plan_dft_1d ( int(n, kind=C_SIZE_T), u, u, &
               MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE+FFTW_MPI_SCRAMBLED_OUT )

        call save_fftw_wisdom

        !call  get_eigenvalues(lambda, 2.0_prec*pi, n, nmin, nmax, PERIODIC)
        !!!!!!!!!!!!!!!!!!!!!!
            m = min(n/2, nmax)
            if (m>=nmin) then
                lambda(nmin:m) = (/ (k, k=nmin-1, m-1) /)
            end if    
            m = max(n/2+1, nmin)
            if (m<=nmax) then
                lambda(m:nmax) = (/ (k, k=m-n-1, nmax-n-1) /)
            end if    
        !!!!!!!!!!!!!!!!!!!!!!
        u = lambda
        call fftw_mpi_execute_dft(plan_backward, u, u)
        call fftw_mpi_execute_dft(plan_forward_scrambled, u, u)
        u = (1.0_prec/real(n,kind=prec))*u
        lambda = real(nint(real(u,kind=prec))**2,kind=prec)
        !lambda = abs(real(u,kind=prec))
        lambda = (-(2.0_prec*pi/d)**2) * lambda

        call fftw_destroy_plan(plan_forward_scrambled)
        call fftw_destroy_plan(plan_backward)
        call fftw_free(c_loc(u))
    end subroutine get_scrambled_eigenvalues
#endif         

#ifdef _QUADPRECISION_
end module tssmq_fourier_common
#else
end module tssm_fourier_common
#endif


