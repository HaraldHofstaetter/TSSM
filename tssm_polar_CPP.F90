
#ifdef _REAL_
 #ifdef _CYLINDRICAL_
  #ifdef _QUADPRECISION_ 
    #define _MODULE_ tssmq_cylindrical_real_3d
  #else
    #define _MODULE_ tssm_cylindrical_real_3d
  #endif  
    #define _METHOD_ cylindrical_real_3d
    #define _WF_ wf_cylindrical_real_3d
 #else
  #ifdef _QUADPRECISION_ 
    #define _MODULE_ tssmq_polar_real_2d
  #else
    #define _MODULE_ tssm_polar_real_2d
  #endif  
    #define _METHOD_ polar_real_2d
    #define _WF_ wf_polar_real_2d
 #endif
 #define _WAVE_FUNCTION_ real_wave_function
 #define _COMPLEX_OR_REAL_ real
#else
 #ifdef _CYLINDRICAL_
  #ifdef _QUADPRECISION_ 
    #define _MODULE_ tssmq_cylindrical_3d
  #else
    #define _MODULE_ tssm_cylindrical_3d
  #endif  
    #define _METHOD_ cylindrical_3d
    #define _WF_ wf_cylindrical_3d
 #else
  #ifdef _QUADPRECISION_ 
    #define _MODULE_ tssmq_polar_2d
  #else
    #define _MODULE_ tssm_polar_2d
  #endif  
    #define _METHOD_ polar_2d
    #define _WF_ wf_polar_2d
 #endif 
 #define _WAVE_FUNCTION_ wave_function
 #define _COMPLEX_OR_REAL_ complex
#endif

#ifdef _QUADPRECISION_ 
#define fftw_alloc_complex fftwq_alloc_complex
#define fftw_alloc_real fftwq_alloc_real
#define fftw_plan_many_dft fftwq_plan_many_dft
#define fftw_plan_many_dft_r2c fftwq_plan_many_dft_r2c
#define fftw_plan_many_dft_c2r fftwq_plan_many_dft_c2r
#define fftw_execute_dft fftwq_execute_dft
#define fftw_execute_dft_r2c fftwq_execute_dft_r2c
#define fftw_execute_dft_c2r fftwq_execute_dft_c2r
#define fftw_destroy_plan fftwq_destroy_plan
#define fftw_free fftwq_free
#endif 


module _MODULE_
#ifdef _QUADPRECISION_ 
    use tssmq_base
    use tssmq_grid
    use tssmq_fourier_common
#else
    use tssm_base
    use tssm_grid
    use tssm_fourier_common
#endif    
    implicit none

    private
    public ::  _METHOD_, _WF_

    type, extends(spectral_method) :: _METHOD_ 
#ifdef _CYLINDRICAL_        
        type(grid_cylindrical_3D) :: g
#else
        type(grid_polar_2D) :: g
#endif
        integer :: nfr
        integer :: nfrmin !integer range for transformed data
        integer :: nfrmax
        integer :: nfthetamin
        integer :: nfthetamax
        integer :: nf1min
        integer :: nf1max
        integer :: nf2min
        integer :: nf2max
#ifdef _CYLINDRICAL_        
        integer :: nfzmin
        integer :: nfzmax
        integer :: nf3min
        integer :: nf3max
#endif        
        real(kind=prec), allocatable :: eigenvalues_r_theta(:,:)
        real(kind=prec), allocatable :: eigenvalues_r(:)
        real(kind=prec), allocatable :: eigenvalues_theta(:)
#ifdef _CYLINDRICAL_        
        real(kind=prec), allocatable :: eigenvalues_z(:)
#endif        
        real(kind=prec), allocatable :: L(:,:,:)
#ifdef _CYLINDRICAL_        
        real(kind=prec), allocatable :: H_z(:,:)
#endif 
        logical :: symmetric_coefficients
#ifdef _OPENMP
        integer, allocatable :: jf(:)
#endif

#ifdef _REAL_
        character(len=32) :: dset_name = "psi_real"  ! name of dataset for hdf5 load/save
#else
        character(len=32) :: dset_name_real = "psi_real" 
        character(len=32) :: dset_name_imag = "psi_imag"
#endif
    contains    
        !final :: final_laguerre_1D
        !! Fortran 2003 feature final seems to be not properly implemented
        !! in the gcc/gfortran compiler :
        procedure :: save => save_method
        procedure :: finalize => finalize_method

    end type _METHOD_

    interface _METHOD_ ! constructor
        module procedure new_method
    end interface _METHOD_


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type, extends(_WAVE_FUNCTION_) ::  _WF_ 
        class(_METHOD_), pointer :: m 

#ifdef _REAL_
        real(kind=prec), pointer  :: up(:)
        complex(kind=prec), pointer  :: ucp(:)
#else
        complex(kind=prec), pointer  :: up(:)
        real(kind=prec), pointer  :: urp(:)
#endif
#ifdef _CYLINDRICAL_        
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf(:,:,:)
#else
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf(:,:)
#endif

#ifdef _REAL_
#ifdef _CYLINDRICAL_
        complex(kind=prec), pointer :: ufc(:,:,:)
#else
        complex(kind=prec), pointer :: ufc(:,:)
#endif
#endif

        type(c_ptr) :: plan_forward 
        type(c_ptr) :: plan_backward
#ifdef _MPI_
        type(c_ptr) :: plan_transpose_forward 
        type(c_ptr) :: plan_transpoes_backward
#endif

        _COMPLEX_OR_REAL_(kind=prec) :: coefficient = 1.0_prec

    contains
        procedure :: create_plans
        procedure :: to_real_space
        procedure :: to_frequency_space
        procedure :: propagate_A
        procedure :: add_apply_A
        procedure :: save
        procedure :: load
        procedure :: set
        procedure :: set_t
#ifndef _REAL_
        procedure :: rset
        procedure :: rset_t
#endif
        procedure :: clone
        procedure :: finalize => finalize_wf
        procedure :: copy
        procedure :: norm
        procedure :: norm_in_frequency_space
        procedure :: normalize
        procedure :: inner_product
        !procedure :: ip_in_frequency_space_ !TODO
        procedure :: distance
        procedure :: axpy
        procedure :: scale
    end type _WF_ 

    interface _WF_ ! constructor
        module procedure new_wf
    end interface _WF_
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains


#ifdef _CYLINDRICAL_
    function new_method(M, nr, nfr, nz, symmetric_coefficients, separated_eigenvalues) result(this)
#else
    function new_method(M, nr, nfr, symmetric_coefficients, separated_eigenvalues) result(this)
#endif
        type(_METHOD_) :: this
        integer, intent(in) :: M
        integer, intent(in) :: nr 
        integer, intent(in) :: nfr 
#ifdef _CYLINDRICAL_
        integer, intent(in) :: nz
#endif
        logical, intent(in), optional :: symmetric_coefficients 
        logical, intent(in), optional :: separated_eigenvalues
        real(kind=prec), parameter :: pi = 3.1415926535897932384626433832795028841971693993751_prec

        integer :: j, j1

!adjust grid parameters, TODO MPI!!!
        this%g%nr = nr ! K+M/2+1
        this%g%nn1min = 0
        this%g%nn1max = nr - 1 ! K+M/2

        this%g%n1min = this%g%nn1min
        this%g%n1max = this%g%nn1max
        this%g%nrmin = this%g%nn1min
        this%g%nrmax = this%g%nn1max
        this%nfr = nfr ! K+1
        this%nfrmin = 0
        this%nfrmax = nfr -1 ! K
        this%nf1min = this%nfrmin
        this%nf1max = this%nfrmax
        this%g%m1min = this%g%nn1min
        this%g%m1max = this%g%nn1max

#ifdef _CYLINDRICAL_
        this%g%nz = nz
        this%g%nn2min = 0
        this%g%nn2max = this%g%nz ! TODO: check this inconsistency !!

        this%g%n2min = this%g%nn2min
        this%g%n2max = this%g%nn2max
        this%g%nzmin = this%g%nn2min
        this%g%nzmax = this%g%nn2max
        this%nf2min = this%g%nn2min
        this%nf2max = this%g%nn2max
        this%nfzmin = this%g%nn2min
        this%nfzmax = this%g%nn2max
        this%g%m2min = this%g%nn2min
        this%g%m2max = this%g%nn2max

        this%g%ntheta = M
        this%g%nn3min = 0
        this%g%nn3max = M - 1

        this%g%n3min = this%g%nn3min
        this%g%n3max = this%g%nn3max
        this%g%nthetamin = this%g%nn3min
        this%g%nthetamax = this%g%nn3max
#ifdef _REAL_
        this%nfthetamin =  0
        this%nfthetamax = this%g%ntheta/2
        this%g%m3min = 0
        this%g%m3max = 2*(this%g%ntheta/2+1)-1
#else
        this%nfthetamin = this%g%nn3min
        this%nfthetamax = this%g%nn3max
        this%g%m3min = this%g%nn3min
        this%g%m3max = this%g%nn3max
#endif
        this%nf3min = this%nfthetamin
        this%nf3max = this%nfthetamax
#else
        this%g%ntheta = M
        this%g%nn2min = 0
        this%g%nn2max = M - 1

        this%g%n2min = this%g%nn2min
        this%g%n2max = this%g%nn2max
        this%g%nthetamin = this%g%nn2min
        this%g%nthetamax = this%g%nn2max
#ifdef _REAL_
        this%nfthetamin =  0
        this%nfthetamax = this%g%ntheta/2
        this%g%m2min = 0
        this%g%m2max = 2*(this%g%ntheta/2+1)-1
#else
        this%nfthetamin = this%g%nn2min
        this%nfthetamax = this%g%nn2max
        this%g%m2min = this%g%nn2min
        this%g%m2max = this%g%nn2max
#endif
        this%nf2min = this%nfthetamin
        this%nf2max = this%nfthetamax
#endif

#ifdef _CYLINDRICAL_
        this%g%alloc_size = max(this%g%m1max-this%g%m1min+1, this%nf1max-this%nf1min+1) &
                            *(this%g%m2max-this%g%m2min+1)*(this%g%m3max-this%g%m3min+1)
#else
        this%g%alloc_size = max(this%g%m1max-this%g%m1min+1, this%nf1max-this%nf1min+1) &
                            *(this%g%m2max-this%g%m2min+1)
#endif

        allocate( this%g%nodes_r(this%g%nrmin:this%g%nrmax ) ) 
        allocate( this%g%weights_r(this%g%nrmin:this%g%nrmax ) ) 
        if (present(symmetric_coefficients).and.symmetric_coefficients) then
            this%symmetric_coefficients = .true.
            allocate( this%L(this%g%nr, this%nfr, 0:M/2 ) )   !TODO: CHECK!!!
        else
            this%symmetric_coefficients = .false.
            allocate( this%L(this%g%nr, this%nfr, 0:M-1 ) )   !TODO: CHECK!!!
        end if

        this%g%dtheta = 2.0_prec*pi/real(M, kind=prec) 
        !this%g%dtheta = 1.0_prec/real(M, kind=prec)  !TODO: CHECK and UNDERSTAND !!!

        allocate( this%g%nodes_theta(this%g%nthetamin:this%g%nthetamax ) ) 
        this%g%nodes_theta = (/ (real(j, kind=prec)*(2.0_prec*pi/real(M, kind=prec)), &
                              j=this%g%nthetamin, this%g%nthetamax ) /)
        
        if ((.not.present(separated_eigenvalues)).or.separated_eigenvalues) then
            allocate( this%eigenvalues_r(this%nfrmin:this%nfrmax ) ) 
            allocate( this%eigenvalues_theta(this%nfthetamin:this%nfthetamax) ) 
        else
            allocate( this%eigenvalues_r_theta(this%nfrmin:this%nfrmax, &
                      this%nfthetamin:this%nfthetamax) ) 
        end if
#ifdef _CYLINDRICAL_
        allocate( this%g%nodes_z(this%g%nzmin:this%g%nzmax ) ) 
        allocate( this%g%weights_z(this%g%nzmin:this%g%nzmax ) ) 
        allocate( this%eigenvalues_z(this%nfzmin:this%nfzmax ) ) 
        allocate( this%H_z(this%g%nzmin:this%g%nzmax, this%nfzmin:this%nfzmax ) ) 
#endif

#ifdef _OPENMP
        allocate( this%g%jj(0:n_threads) )
        allocate( this%jf(0:n_threads) )
        do j=0,n_threads-1
           this%g%jj(j) = j*ceiling(real(this%g%nthetamax-this%g%nthetamin+1)/real(n_threads))
           this%jf(j) = j*ceiling(real(this%nfthetamax-this%nfthetamin+1)/real(n_threads))
        end do
        this%g%jj(n_threads) = this%g%nthetamax-this%g%nthetamin+1
        this%jf(n_threads) = this%nfthetamax-this%nfthetamin+1 
#endif

    end function new_method


    subroutine finalize_method(this)
        class(_METHOD_), intent(inout) :: this

        deallocate( this%g%nodes_r )
        deallocate( this%g%weights_r )
        deallocate( this%L )
        if (allocated(this%eigenvalues_r_theta)) then
            deallocate( this%eigenvalues_r_theta )
        else
            deallocate( this%eigenvalues_r )
            deallocate( this%eigenvalues_theta )
        endif
#ifdef _CYLINDRICAL_
        deallocate( this%g%nodes_z )
        deallocate( this%g%weights_z )
        deallocate( this%H_z )
        deallocate( this%eigenvalues_z )
#endif
#ifdef _OPENMP
        deallocate( this%g%jj )
        deallocate( this%jf )
#endif        
    end subroutine finalize_method

    subroutine save_method(this, filename)
#ifdef _QUADPRECISION_ 
        use tssmq_hdf5_helper
#else
        use tssm_hdf5_helper
#endif
        class(_METHOD_), intent(inout) :: this
        character(len=*), intent(in) :: filename
           
        call create_file(filename)
        call write_integer_attr(filename, "nr",  this%g%nr)
        call write_integer_attr(filename, "nfr",  this%nfr)
        call write_integer_attr(filename, "ntheta",  this%g%ntheta)
        call write_array(filename, "nodes_r", this%g%nodes_r, 1, (/ this%g%nr /) )
        call write_array(filename, "weights_r", this%g%weights_r, 1, (/ this%g%nr /) )
        if (allocated(this%eigenvalues_r_theta)) then
            call write_array(filename, "eigenvalues_r_theta", this%eigenvalues_r_theta, &
                     2, (/ this%nfr, this%g%ntheta/2+1 /) )
        else
            call write_array(filename, "eigenvalues_r", this%eigenvalues_r, &
                     1, (/ this%nfr /) )
            call write_array(filename, "eigenvalues_theta", this%eigenvalues_theta, &
                     1, (/ this%g%ntheta /) )
        endif
        call write_array(filename, "L", this%L, 3, (/ this%g%nr, this%nfr, this%g%ntheta/2+1 /) )

#ifdef _CYLINDRICAL_
        call write_integer_attr(filename, "nz",  this%g%nz)
        call write_array(filename, "nodes_z", this%g%nodes_z, 1, (/ this%g%nz /) )
        call write_array(filename, "weights_z", this%g%weights_z, 1, (/ this%g%nz /) )
        call write_array(filename, "eigenvalues_z", this%eigenvalues_z, &
                     1, (/ this%g%nz /) )
        call write_array(filename, "H_z", this%H_z, 2, (/ this%g%nz, this%g%nz /) )
#endif
    end subroutine save_method


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function new_wf(m, u, coefficient) result(this)
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_ptr, c_loc
        type(_WF_) :: this
        class(_METHOD_), target, intent(inout) :: m
#ifdef _REAL_       
#ifndef _CYLINDRICAL_
        real(kind=prec), optional, target, intent(inout) :: u(:,:)
#else
        real(kind=prec), optional, target, intent(inout) :: u(:,:,:)
#endif
        real(kind=prec), optional, intent(in) :: coefficient
#else            
#ifndef _CYLINDRICAL_
        complex(kind=prec), optional, target, intent(inout) :: u(:,:)
#else            
        complex(kind=prec), optional, target, intent(inout) :: u(:,:,:)
#endif
        complex(kind=prec), optional, intent(in) :: coefficient
#endif

#ifdef _REAL_       
#ifndef _CYLINDRICAL_
        real(kind=prec), pointer  :: umem(:,:)
#else            
        real(kind=prec), pointer  :: umem(:,:,:)
#endif
#else            
#ifndef _CYLINDRICAL_
        complex(kind=prec), pointer  :: umem(:,:)
#else            
        complex(kind=prec), pointer  :: umem(:,:,:)
#endif
#endif

        type(C_PTR) :: p

        call initialize_tssm_fourier

        this%m => m
        if (present(coefficient)) then
            this%coefficient = coefficient
        end if    

        if (present(u)) then
#ifndef _CYLINDRICAL_
            umem => u(m%g%m1min:m%g%m1max, m%g%m2min:m%g%m2max) 
#else
            umem => u(m%g%m1min:m%g%m1max, m%g%m2min:m%g%m2max, m%g%m3min:m%g%m3max) 
#endif
        else
#ifdef _REAL_        
#ifdef _QUADPRECISION_
            p = fftw_alloc_complex(this%m%g%alloc_size/2)
#else
            p = fftw_alloc_real(this%m%g%alloc_size)
#endif 
#else            
            p = fftw_alloc_complex(this%m%g%alloc_size)
#endif
            call c_f_pointer(p, this%up, [this%m%g%alloc_size])
#ifdef _REAL_        
            call c_f_pointer(p, this%ucp, [this%m%g%alloc_size/2])
#else
            call c_f_pointer(p, this%urp, [2*this%m%g%alloc_size])
#endif 
#ifndef _CYLINDRICAL_
            umem(this%m%g%m1min:this%m%g%m1max, this%m%g%m2min:this%m%g%m2max) => this%up
            this%u => umem(this%m%g%nrmin:this%m%g%nrmax, this%m%g%nthetamin:this%m%g%nthetamax)
#else
            umem(this%m%g%m1min:this%m%g%m1max, this%m%g%m2min:this%m%g%m2max, this%m%g%m3min:this%m%g%m3max) => this%up 
            this%u => umem(this%m%g%nrmin:this%m%g%nrmax, this%m%g%nzmin:this%m%g%nzmax, this%m%g%nthetamin:this%m%g%nthetamax)  
#endif
        end if

#ifndef _CYLINDRICAL_
        this%uf(this%m%g%nrmin:this%m%g%nrmax, this%m%nfthetamin:this%m%nfthetamax) => this%up
#else
        this%uf(this%m%g%nrmin:this%m%g%nrmax, this%m%nfzmin:this%m%nfzmax, this%m%nfthetamin:this%m%nfthetamax) =>this%up  
#endif
         call this%create_plans
#ifdef _REAL_        
#ifndef _CYLINDRICAL_
        this%ufc(this%m%g%nrmin:this%m%g%nrmax, this%m%nfthetamin:this%m%nfthetamax) => this%ucp
#else
        this%ufc(this%m%g%nrmin:this%m%g%nrmax, this%m%nfzmin:this%m%nfzmax, this%m%nfthetamin:this%m%nfthetamax) =>this%ucp  
#endif
#endif
    end function new_wf


    
    subroutine create_plans(this)
        class(_WF_) :: this
#ifdef _MPI_
        integer(kind=C_SIZE_T) :: dft_dim(1)  
        integer(kind=C_SIZE_T) :: howmany
#else
        integer :: dft_dim(1)  
        integer :: howmany
#endif
        

        call load_fftw_wisdom
#ifdef _MPI_
     !       this%plan_forward_transpose = fftw_mpi_plan_transpose(int(this%m%g%ntheta, kind=C_SIZE_T), &
     !            int(this%m%g%nr, kind=C_SIZE_T), this%up, this%up, MPI_COMM_WORLD, &
     !            fftw_planning_rigor + FFTW_MPI_TRANSPOSED_OUT)
     !       this%plan_backward_transpose = fftw_mpi_plan_transpose(int(this%m%g%ntheta, kind=C_SIZE_T), &
     !            int(this%m%g%nr, kind=C_SIZE_T), this%up, this%up, MPI_COMM_WORLD, &
     !            fftw_planning_rigor + FFTW_MPI_TRANSPOSED_IN)

#else
        dft_dim = (/ this%m%g%ntheta /)
#ifdef _CYLINDRICAL_
        howmany = this%m%g%nr*this%m%g%nz
#else        
        howmany = this%m%g%nr
#endif        


#ifdef _REAL_
        this%plan_forward = fftw_plan_many_dft_r2c(1, dft_dim, howmany, &   ! rank, *n, howmany,
                                               this%up, &                   ! *in, 
                                               dft_dim, this%m%g%nr, 1, &   ! *inembed, istride, idist,
                                               this%ucp, &                  ! *out, 
                                               dft_dim, this%m%g%nr, 1, &   ! *onembed, ostride, odist,
                                               fftw_planning_rigor)         ! flags
        this%plan_backward = fftw_plan_many_dft_c2r(1, dft_dim, howmany, &  ! rank, *n, howmany,
                                               this%ucp, &                  ! *in, 
                                               dft_dim, this%m%g%nr, 1, &   ! *inembed, istride, idist,
                                               this%up, &                   ! *out, 
                                               dft_dim, this%m%g%nr, 1, &   ! *onembed, ostride, odist,
                                               fftw_planning_rigor)         ! flags
#else
        this%plan_forward = fftw_plan_many_dft(1, dft_dim, howmany, &       ! rank, *n, howmany,
                                               this%up, &                   ! *in, 
                                               dft_dim, this%m%g%nr, 1, &   ! *inembed, istride, idist,
                                               this%up, &                   ! *out, 
                                               dft_dim, this%m%g%nr, 1, &   ! *onembed, ostride, odist,
                                               FFTW_FORWARD, fftw_planning_rigor) ! sign, flags
        this%plan_backward = fftw_plan_many_dft(1, dft_dim, howmany, &      ! rank, *n, howmany,
                                               this%up, &                   ! *in, 
                                               dft_dim, this%m%g%nr, 1, &   ! *inembed, istride, idist,
                                               this%up, &                   ! *out, 
                                               dft_dim, this%m%g%nr, 1, &   ! *onembed, ostride, odist,
                                               FFTW_BACKWARD, fftw_planning_rigor) ! sign, flags
#endif            

#endif            
        call save_fftw_wisdom
    end subroutine create_plans



    subroutine to_real_space(this)
        class(_WF_), intent(inout) :: this
        integer :: m, m1
        real(kind=prec), parameter :: pi2 = 2.0_prec*3.1415926535897932384626433832795028841971693993751_prec

        if (.not. this%is_real_space) then
!$OMP PARALLEL DO PRIVATE(m, m1)
            do m = this%m%nfthetamin, this%m%nfthetamax
                m1 = m
                if (this%m%symmetric_coefficients .and. (m1>this%m%g%ntheta/2)) then
                    m1 = abs(m1-this%m%g%ntheta)
                end if
#ifdef _REAL_
#ifndef _CYLINDRICAL_
                this%ufc(:,m) = matmul((this%m%L(:,:,m1)), this%ufc(this%m%nfrmin:this%m%nfrmax,m))
#else
                this%ufc(this%m%nfrmin:this%m%nfrmax,:,m) = matmul(this%ufc(this%m%nfrmin:this%m%nfrmax,:,m), &
                                                                  transpose(this%m%H_z))
                this%ufc(:,:,m) = matmul((this%m%L(:,:,m1)), this%ufc(this%m%nfrmin:this%m%nfrmax,:,m))
#endif
#else
#ifndef _CYLINDRICAL_
                this%uf(:,m) = matmul((this%m%L(:,:,m1)), this%uf(this%m%nfrmin:this%m%nfrmax,m))
#else
                this%uf(this%m%nfrmin:this%m%nfrmax,:,m) = matmul(this%uf(this%m%nfrmin:this%m%nfrmax,:,m), &
                                                                  transpose(this%m%H_z))
                this%uf(:,:,m) = matmul((this%m%L(:,:,m1)), this%uf(this%m%nfrmin:this%m%nfrmax,:,m))
#endif
#endif
            end do    
!$OMP END PARALLEL DO            
#ifdef _MPI_
#ifdef _REAL_
        !TODO
!            call fftw_mpi_execute_dft_r2c(this%plan_forward, this%up, this%ucp)
#else
        !TODO
!            call fftw_mpi_execute_dft(this%plan_forward, this%up, this%up)
#endif
#else
#ifdef _REAL_
            call fftw_execute_dft_c2r(this%plan_backward, this%ucp, this%up)
#else
            call fftw_execute_dft(this%plan_backward, this%up, this%up)
#endif
#endif
           this%is_real_space = .true.
        end if
    end subroutine to_real_space


    subroutine to_frequency_space(this)
        class(_WF_), intent(inout) :: this
        integer :: m, m1
        real(kind=prec) :: f
        real(kind=prec), parameter :: pi2 = 2.0_prec*3.1415926535897932384626433832795028841971693993751_prec

        if (this%is_real_space) then
#ifdef _MPI_
#ifdef _REAL_
        !TODO
!            call fftw_mpi_execute_dft_r2c(this%plan_forward, this%up, this%ucp)
#else
        !TODO
!            call fftw_mpi_execute_dft(this%plan_forward, this%up, this%up)
#endif
#else
#ifdef _REAL_
            call fftw_execute_dft_r2c(this%plan_forward, this%up, this%ucp)
#else
            call fftw_execute_dft(this%plan_forward, this%up, this%up)
#endif
#endif
            f = pi2/this%m%g%ntheta
!$OMP PARALLEL DO PRIVATE(m, m1)
            do m = this%m%nfthetamin, this%m%nfthetamax
                m1 = m
                if (this%m%symmetric_coefficients .and. (m1>this%m%g%ntheta/2)) then
                    m1 = abs(m1-this%m%g%ntheta)
                end if
#ifdef _REAL_
#ifndef _CYLINDRICAL_
                this%ufc(this%m%nfrmin:this%m%nfrmax,m) = matmul(transpose(this%m%L(:,:,m1)), &
                                                         this%m%g%weights_r*this%ufc(:,m))
                this%ufc(this%m%nfrmin:this%m%nfrmax,m) = f*this%ufc(this%m%nfrmin:this%m%nfrmax,m)
#else
                this%ufc(this%m%nfrmin:this%m%nfrmax,:,m) = matmul( transpose(this%m%L(:,:,m1)), &
                              spread(this%m%g%weights_r, 2, this%m%nfzmax-this%m%nfzmin+1) * this%ufc(:,:,m) )
                this%ufc(this%m%nfrmin:this%m%nfrmax,:,m) = f*this%ufc(this%m%nfrmin:this%m%nfrmax,:,m)
                this%ufc(this%m%nfrmin:this%m%nfrmax,:,m) = matmul(spread(this%m%g%weights_z, 1, &
                              this%m%nfrmax-this%m%nfrmin+1)*this%ufc(this%m%nfrmin:this%m%nfrmax,:,m), this%m%H_z)
#endif

#else
#ifndef _CYLINDRICAL_
                this%uf(this%m%nfrmin:this%m%nfrmax,m) = matmul(transpose(this%m%L(:,:,m1)), &
                                                         this%m%g%weights_r*this%uf(:,m))
                this%uf(this%m%nfrmin:this%m%nfrmax,m) = f*this%uf(this%m%nfrmin:this%m%nfrmax,m)
#else
                this%uf(this%m%nfrmin:this%m%nfrmax,:,m) = matmul( transpose(this%m%L(:,:,m1)), &
                              spread(this%m%g%weights_r, 2, this%m%nfzmax-this%m%nfzmin+1) * this%uf(:,:,m) )
                this%uf(this%m%nfrmin:this%m%nfrmax,:,m) = f*this%uf(this%m%nfrmin:this%m%nfrmax,:,m)
                this%uf(this%m%nfrmin:this%m%nfrmax,:,m) = matmul(spread(this%m%g%weights_z, 1, &
                              this%m%nfrmax-this%m%nfrmin+1)*this%uf(this%m%nfrmin:this%m%nfrmax,:,m), this%m%H_z)
#endif
#endif
            end do    
!$OMP END PARALLEL DO            
           this%is_real_space = .false.
        end if    
    end subroutine to_frequency_space
   

    subroutine propagate_A(this, dt)
        class(_WF_), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt
#ifndef _CYLINDRICAL_
        complex(kind=prec), pointer :: uf(:,:)
#else
        complex(kind=prec), pointer :: uf(:,:,:)
#endif 
#ifdef _OPENMP
        integer :: j
        real(kind=prec), pointer :: ev(:)
        real(kind=prec), pointer :: ev1(:,:)
#endif        

        call this%to_frequency_space

#ifndef _CYLINDRICAL_
#ifndef _OPENMP
#ifdef _REAL_
        uf => this%ufc(this%m%nfrmin:this%m%nfrmax,:)
#else
        uf => this%uf(this%m%nfrmin:this%m%nfrmax,:)
#endif        
        if (allocated(this%m%eigenvalues_r_theta)) then
            uf = exp((dt*this%coefficient)*this%m%eigenvalues_r_theta) * uf
        else 
            uf = spread(exp((dt*this%coefficient)*this%m%eigenvalues_r), &
                        2, this%m%nfthetamax-this%m%nfthetamin+1) * uf
            uf = spread(exp((dt*this%coefficient)*this%m%eigenvalues_theta), &
                        1, this%m%nfrmax-this%m%nfrmin+1) * uf
        end if
#else 
!$OMP PARALLEL DO PRIVATE(j, uf, ev, ev1) 
        do j=1,n_threads
#ifdef _REAL_
            uf => this%ufc(this%m%nfrmin:this%m%nfrmax, &
                           lbound(this%uf,2)+this%m%jf(j-1):&
                           lbound(this%uf,2)+this%m%jf(j)-1)
#else            
            uf => this%uf(this%m%nfrmin:this%m%nfrmax, &
                          lbound(this%uf,2)+this%m%jf(j-1):&
                          lbound(this%uf,2)+this%m%jf(j)-1)
#endif            
          if (allocated(this%m%eigenvalues_r_theta)) then
            ev1 => this%m%eigenvalues_r_theta(:, &
                          lbound(this%m%eigenvalues_r_theta,2)+this%m%jf(j-1):&
                          lbound(this%m%eigenvalues_r_theta,2)+this%m%jf(j)-1)
            uf = exp((dt*this%coefficient)*ev1) * uf
          else
            ev => this%m%eigenvalues_theta( &
                          lbound(this%m%eigenvalues_theta,1)+this%m%jf(j-1):&
                          lbound(this%m%eigenvalues_theta,1)+this%m%jf(j)-1)
            uf = spread(exp((dt*this%coefficient)*this%m%eigenvalues_r),&
                          2, this%m%jf(j)-this%m%jf(j-1)) * uf
            uf = spread(exp((dt*this%coefficient)*ev), &
                          1, this%m%nfrmax-this%m%nfrmin+1) * uf
          end if  
        end do
!$OMP END PARALLEL DO 

#endif

#else        

#ifndef _OPENMP
#ifdef _REAL_
        uf => this%ufc(this%m%nfrmin:this%m%nfrmax,:,:)
#else
        uf => this%uf(this%m%nfrmin:this%m%nfrmax,:,:)
#endif
        if (allocated(this%m%eigenvalues_r_theta)) then
          uf = spread(exp((dt*this%coefficient)*this%m%eigenvalues_r_theta), &
                      2, this%m%nfzmax-this%m%nfzmin+1 )
        else
          uf = spread(spread(exp((dt*this%coefficient)*this%m%eigenvalues_r), &
                      2, this%m%nfzmax-this%m%nfzmin+1), &
                      3, this%m%nfthetamax-this%m%nfthetamin+1) * uf
          uf = spread(spread(exp((dt*this%coefficient)*this%m%eigenvalues_theta), &
                      1, this%m%nfrmax-this%m%nfrmin+1), & 
                      2, this%m%nfzmax-this%m%nfzmin+1) * uf
        endif                         
        uf = spread(spread(exp((dt*this%coefficient)*this%m%eigenvalues_z), &
                    1, this%m%nfrmax-this%m%nfrmin+1), & 
                    3, this%m%nfthetamax-this%m%nfthetamin+1) * uf
#else 
!$OMP PARALLEL DO PRIVATE(j, uf, ev, ev1) 
        do j=1,n_threads
#ifdef _REAL_
            uf => this%ufc(this%m%nfrmin:this%m%nfrmax,:,&
                           lbound(this%uf,3)+this%m%jf(j-1):&
                           lbound(this%uf,3)+this%m%jf(j)-1)
#else
            uf => this%uf(this%m%nfrmin:this%m%nfrmax,:,&
                          lbound(this%uf,3)+this%m%jf(j-1):&
                          lbound(this%uf,3)+this%m%jf(j)-1)
#endif            
            if (allocated(this%m%eigenvalues_r_theta)) then
                ev1 => this%m%eigenvalues_r_theta(:, &
                           lbound(this%m%eigenvalues_r_theta,1)+this%m%jf(j-1):&
                           lbound(this%m%eigenvalues_r_theta,1)+this%m%jf(j)-1)
                uf = spread(exp((dt*this%coefficient)*ev1), &
                           2, this%m%nfzmax-this%m%nfzmin+1 )
            else
                ev => this%m%eigenvalues_theta(&
                           lbound(this%m%eigenvalues_theta,1)+this%m%jf(j-1):&
                           lbound(this%m%eigenvalues_theta,1)+this%m%jf(j)-1)
                uf = spread(spread(exp((dt*this%coefficient)*this%m%eigenvalues_r), &
                           2, this%m%nfzmax-this%m%nfzmin+1), &
                           3,  this%m%jf(j)-this%m%jf(j-1)) * uf
                uf = spread(spread(exp((dt*this%coefficient)*ev), &
                           1, this%m%nfrmax-this%m%nfrmin+1), &
                           2, this%m%nfzmax-this%m%nfzmin+1) * uf
            endif                
            uf = spread(spread(exp((dt*this%coefficient)*this%m%eigenvalues_z), &
                           1, this%m%nfrmax-this%m%nfrmin+1), &
                           3,  this%m%jf(j)-this%m%jf(j-1)) * uf
        end do
!$OMP END PARALLEL DO 
#endif
#endif

    end subroutine propagate_A


    subroutine add_apply_A(this, wf, coefficient)
        class(_WF_), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: wf 
        _COMPLEX_OR_REAL_(kind=prec), intent(in), optional :: coefficient
       _COMPLEX_OR_REAL_(kind=prec) :: C 
#ifndef _CYLINDRICAL_
        complex(kind=prec), pointer :: uf1(:,:)
        complex(kind=prec), pointer :: uf2(:,:)
#else
        complex(kind=prec), pointer :: uf1(:,:,:)
        complex(kind=prec), pointer :: uf2(:,:,:)
#endif 
#ifdef _OPENMP
        integer :: j
        real(kind=prec), pointer :: ev(:)
        real(kind=prec), pointer :: ev1(:,:)
#endif        

        select type (wf)
        class is (_WF_)
        
        if (.not.associated(wf%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    

        C = this%coefficient
        if (present(coefficient)) then
            C = C*coefficient
        end if    

        call this%to_frequency_space
        call wf%to_frequency_space

#ifndef _CYLINDRICAL_
#ifndef _OPENMP
#ifdef _REAL_
        uf1 => this%ufc(this%m%nfrmin:this%m%nfrmax,:)
        uf2 => wf%ufc(this%m%nfrmin:this%m%nfrmax,:)
#else
        uf1 => this%uf(this%m%nfrmin:this%m%nfrmax,:)
        uf2 => wf%uf(this%m%nfrmin:this%m%nfrmax,:)
#endif        
        if (allocated(this%m%eigenvalues_r_theta)) then
            uf2 = uf2 + (C*this%m%eigenvalues_r_theta) * uf1
        else
            uf2 = uf2 + spread(C*this%m%eigenvalues_r, &
                      2, this%m%nfthetamax-this%m%nfthetamin+1) * uf1
            uf2 = uf2 + spread(C*this%m%eigenvalues_theta, &
                      1, this%m%nfrmax-this%m%nfrmin+1) * uf1
        endif    
#else 
!$OMP PARALLEL DO PRIVATE(j, uf1, uf2, ev, ev1) 
        do j=1,n_threads
#ifdef _REAL_
          uf1 => this%ufc(this%m%nfrmin:this%m%nfrmax, &
                 lbound(this%uf,2)+this%m%jf(j-1):&
                 lbound(this%uf,2)+this%m%jf(j)-1)
          uf2 => wf%ufc(this%m%nfrmin:this%m%nfrmax,&
                 lbound(this%uf,2)+this%m%jf(j-1):&
                 lbound(this%uf,2)+this%m%jf(j)-1)
#else            
          uf1 => this%uf(this%m%nfrmin:this%m%nfrmax,&
                 lbound(this%uf,2)+this%m%jf(j-1):&
                 lbound(this%uf,2)+this%m%jf(j)-1)
          uf2 => wf%uf(this%m%nfrmin:this%m%nfrmax,&
                 lbound(this%uf,2)+this%m%jf(j-1):&
                 lbound(this%uf,2)+this%m%jf(j)-1)
#endif            
          if (allocated(this%m%eigenvalues_r_theta)) then
            ev1 => this%m%eigenvalues_r_theta(:, &
                 lbound(this%m%eigenvalues_r_theta,2)+this%m%jf(j-1):&
                 lbound(this%m%eigenvalues_r_theta,2)+this%m%jf(j)-1)
            uf2 = uf2 + (C*ev1) * uf1
          else
            ev => this%m%eigenvalues_theta( &
                 lbound(this%m%eigenvalues_theta,1)+this%m%jf(j-1):&
                 lbound(this%m%eigenvalues_theta,1)+this%m%jf(j)-1)
            uf2 = uf2 + spread(C*this%m%eigenvalues_r, &
                 2, this%m%jf(j)-this%m%jf(j-1)) * uf1
            uf2 = uf2 + spread(C*ev ,&
                 1, this%m%nfrmax-this%m%nfrmin+1) * uf1
          endif  
        end do
!$OMP END PARALLEL DO 

#endif

#else        

#ifndef _OPENMP
#ifdef _REAL_
        uf1 => this%ufc(this%m%nfrmin:this%m%nfrmax,:,:)
        uf2 => wf%ufc(this%m%nfrmin:this%m%nfrmax,:,:)
#else
        uf1 => this%uf(this%m%nfrmin:this%m%nfrmax,:,:)
        uf2 => wf%uf(this%m%nfrmin:this%m%nfrmax,:,:)
#endif
        if (allocated(this%m%eigenvalues_r_theta)) then
        else
          uf2 = uf2 + spread(spread(C*this%m%eigenvalues_r, &
                  2, this%m%nfzmax-this%m%nfzmin+1), &
                  3, this%m%nfthetamax-this%m%nfthetamin+1) * uf1
          uf2 = uf2 + spread(spread(C*this%m%eigenvalues_theta, &
                  1, this%m%nfrmax-this%m%nfrmin+1), &
                  2, this%m%nfzmax-this%m%nfzmin+1) * uf1
        endif          
        uf2 = uf2 + spread(spread(C*this%m%eigenvalues_z, &
                  1, this%m%nfrmax-this%m%nfrmin+1), &
                  3, this%m%nfthetamax-this%m%nfthetamin+1) * uf1
#else 
!$OMP PARALLEL DO PRIVATE(j, uf1, uf2, ev, ev1) 
        do j=1,n_threads
#ifdef _REAL_
            uf1 => this%ufc(this%m%nfrmin:this%m%nfrmax,:,&
                  lbound(this%uf,3)+this%m%jf(j-1):&
                  lbound(this%uf,3)+this%m%jf(j)-1)
            uf2 => wf%ufc(this%m%nfrmin:this%m%nfrmax,:,&
                  lbound(this%uf,3)+this%m%jf(j-1):&
                  lbound(this%uf,3)+this%m%jf(j)-1)
#else
            uf1 => this%uf(this%m%nfrmin:this%m%nfrmax,:,&
                  lbound(this%uf,3)+this%m%jf(j-1):&
                  lbound(this%uf,3)+this%m%jf(j)-1)
            uf2 => wf%uf(this%m%nfrmin:this%m%nfrmax,:,&
                  lbound(this%uf,3)+this%m%jf(j-1):&
                  lbound(this%uf,3)+this%m%jf(j)-1)
#endif            
            if (allocated(this%m%eigenvalues_r_theta)) then
              ev1 => this%m%eigenvalues_r_theta(:, &
                  lbound(this%m%eigenvalues_r_theta,2)+this%m%jf(j-1):&
                  lbound(this%m%eigenvalues_r_theta,2)+this%m%jf(j)-1)
              uf2 = uf2 + spread(C*ev1, 2, this%m%nfzmax-this%m%nfzmin+1 ) * uf1
            else
              ev => this%m%eigenvalues_theta(&
                  lbound(this%m%eigenvalues_theta,1)+this%m%jf(j-1):&
                  lbound(this%m%eigenvalues_theta,1)+this%m%jf(j)-1)
              uf2 = uf2 + spread(spread(C*this%m%eigenvalues_r, &
                  2, this%m%nfzmax-this%m%nfzmin+1), &
                  3,  this%m%jf(j)-this%m%jf(j-1)) * uf1
              uf2 = uf2 + spread(spread(C*ev, &
                  1, this%m%nfrmax-this%m%nfrmin+1), &
                  2, this%m%nfzmax-this%m%nfzmin+1) * uf1
            endif    
            uf2 = uf2 + spread(spread(C*this%m%eigenvalues_z, &
                  1, this%m%nfrmax-this%m%nfrmin+1), &
                  3,  this%m%jf(j)-this%m%jf(j-1)) * uf1
        end do
!$OMP END PARALLEL DO 
#endif
#endif
        class default
           stop "E: wave functions not belonging to the same method"
        end select
  
    end subroutine add_apply_A



    subroutine save(this, filename)
#ifdef _NO_HDF5_
        class(_WF_), intent(inout) :: this
        character(len=*), intent(in) :: filename
        print *, "W: save not implemented"
#else
#ifdef _QUADPRECISION_
        use tssmq_hdf5
#else
        use tssm_hdf5
#endif        
        class(_WF_), intent(inout) :: this
        character(len=*), intent(in) :: filename
        !TODO
       call this%to_real_space

#ifdef _REAL_        
        call hdf5_save_real_gridfun(this%m%g, this%u, filename, trim(this%m%dset_name))
#else        
        call hdf5_save_complex_gridfun(this%m%g, this%u, filename, &
                     trim(this%m%dset_name_real), trim(this%m%dset_name_imag))
#endif        

#ifndef _CYLINDRICAL_
        call hdf5_write_attributes(filename, (/ "grid" /), &
          (/ "polar_2D" /), &
          (/ "nr    ", "ntheta" /), (/ this%m%g%nr, this%m%g%ntheta /), &
          (/ character(len=0) :: /), (/ real(prec) :: /) ) 
#else        
        call hdf5_write_attributes(filename, (/ "grid" /), &
          (/ "cylindrical_3d"  /), &
          (/ "nr    ", "ntheta", "nz    " /), (/ this%m%g%nr, this%m%g%ntheta, this%m%g%nz /), &
          (/ character(len=0) :: /), (/ real(prec) :: /) ) 
#endif        
!TODO save nodes !!!
#endif        
    end subroutine save


   subroutine load(this, filename)
#ifdef _NO_HDF5_
        class(_WF_), intent(inout) :: this
        character(len=*), intent(in) :: filename
        print *, "W: load not implemented"
#else
#ifdef _QUADPRECISION_
        use tssmq_hdf5
#else
        use tssm_hdf5
#endif        
        class(_WF_), intent(inout) :: this
        character(len=*), intent(in) :: filename
#ifdef _QUADPRECISION_
        real(kind=prec), parameter :: eps = epsilon(1.0_8)
#else        
        real(kind=prec), parameter :: eps = epsilon(1.0_prec)
#endif
        character(len=20) :: s
        integer :: n, len

!TODO check compatibility of grids!!!     
       call hdf5_read_string_attribute(filename, "grid", s, len)      
#ifndef _CYLINDRICAL_
       if (s(1:len)/="polar_2D") then
              stop "E: wrong grid type"
       end if
#else
       if (s(1:len)/="cylindrical_3D") then
              stop "E: wrong grid type"
       end if
#endif
       if ((hdf5_read_integer_attribute(filename, "nr") /= this%m%g%nr)) then
              stop "E: incompatible grids"
       end if       
       if (hdf5_read_integer_attribute(filename, "ntheta") /= this%m%g%ntheta) then
              stop "E: incompatible grids"
       end if       
#ifdef _CYLINDRICAL_
       if ((hdf5_read_integer_attribute(filename, "nz") /= this%m%g%nz)) then 
              stop "E: incompatible grids"
       end if       
#endif

#ifdef _REAL_        
        call hdf5_load_real_gridfun(this%m%g, this%u, filename, trim(this%m%dset_name))
#else      
        call hdf5_load_complex_gridfun(this%m%g, this%u, filename, & 
                     trim(this%m%dset_name_real), trim(this%m%dset_name_imag))
#endif        
        this%is_real_space = .true.
#endif        

    end subroutine load


    subroutine set(this, f)
        class(_WF_), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), external :: f
#ifdef _REAL_        
        call this%m%g%set_real_gridfun(this%u, f)
#else
        call this%m%g%set_complex_gridfun(this%u, f)
#endif        
        this%is_real_space = .true.
    end subroutine set


    subroutine set_t(this, f, t)
        class(_WF_), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), external :: f
        real(kind=prec), intent(in) :: t
#ifdef _REAL_        
        call this%m%g%set_t_real_gridfun(this%u, f, t)
#else
        call this%m%g%set_t_complex_gridfun(this%u, f, t)
#endif        
        this%is_real_space = .true.
    end subroutine set_t


#ifndef _REAL_        

    subroutine rset(this, f)
        class(_WF_), intent(inout) :: this
        real(kind=prec), external :: f
        call this%m%g%rset_complex_gridfun(this%u, f)
        this%is_real_space = .true.
    end subroutine rset

    subroutine rset_t(this, f, t)
        class(_WF_), intent(inout) :: this
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t
        call this%m%g%rset_t_complex_gridfun(this%u, f, t)
        this%is_real_space = .true.
    end subroutine rset_t

#endif        


    function clone(this) 
        class(_WF_), intent(inout) :: this
        class(_WAVE_FUNCTION_), pointer :: clone
        type(_WF_), pointer :: p

        allocate( p )
        p = _WF_(this%m, coefficient=this%coefficient)
        clone => p
        
    end function clone


    subroutine finalize_wf(this)
        use, intrinsic :: iso_c_binding, only: c_loc
        class(_WF_), intent(inout) :: this

        call fftw_destroy_plan(this%plan_forward)
        call fftw_destroy_plan(this%plan_backward)
       !call fftw_free(c_loc(this%up))
        call fftw_free(c_loc(this%up(1)))
    end subroutine finalize_wf


    subroutine copy(this, source)
        class(_WF_), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: source 
#ifdef _OPENMP
#ifndef _CYLINDRICAL_
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u1(:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u2(:,:)
#ifdef _REAL_
        complex(kind=prec), pointer :: ufc1(:,:)
        complex(kind=prec), pointer :: ufc2(:,:)
#else
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf1(:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf2(:,:)
#endif
#else
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u1(:,:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u2(:,:,:)
#ifdef _REAL_
        complex(kind=prec), pointer :: ufc1(:,:,:)
        complex(kind=prec), pointer :: ufc2(:,:,:)
#else
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf1(:,:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf2(:,:,:)
#endif
#endif
        integer :: j
#endif        
        select type (source)
        class is (_WF_)
        if (.not.associated(source%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    
        if (this%is_real_space) then
#ifndef _OPENMP
            this%u = source%u
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
           do j=1,n_threads
#ifndef _CYLINDRICAL_
              u1 => this%u(:,lbound(this%u,2)+this%m%g%jj(j-1):lbound(this%u,2)+this%m%g%jj(j)-1)
              u2 => source%u(:,lbound(source%u,2)+this%m%g%jj(j-1):lbound(source%u,2)+this%m%g%jj(j)-1)
#else
              u1 => this%u(:,:,lbound(this%u,3)+this%m%g%jj(j-1):lbound(this%u,3)+this%m%g%jj(j)-1)
              u2 => source%u(:,:,lbound(source%u,3)+source%m%g%jj(j-1):lbound(source%u,3)+this%m%g%jj(j)-1)
#endif
              u1 = u2
           end do
!$OMP END PARALLEL DO 
#endif
        else    
#ifdef _REAL_
#ifndef _OPENMP
          this%ufc = source%ufc
#else
!$OMP PARALLEL DO PRIVATE(j, ufc1, ufc2) 
          do j=1,n_threads
#ifndef _CYLINDRICAL_
              ufc1 => this%ufc(:,lbound(this%ufc,2)+this%m%jf(j-1):lbound(this%ufc,2)+this%m%jf(j)-1)
              ufc2 => source%ufc(:,lbound(source%ufc,2)+this%m%jf(j-1):lbound(source%ufc,2)+this%m%jf(j)-1)
#else
              ufc1 => this%ufc(:,:,lbound(this%ufc,3)+this%m%jf(j-1):lbound(this%ufc,3)+this%m%jf(j)-1)
              ufc2 => source%ufc(:,:,lbound(source%ufc,3)+source%m%jf(j-1):lbound(source%ufc,3)+this%m%jf(j)-1)
#endif
              ufc1 = ufc2
          end do
!$OMP END PARALLEL DO 
#endif

#else
#ifndef _OPENMP
                this%uf = source%uf 
#else
!$OMP PARALLEL DO PRIVATE(j, uf1, uf2) 
             do j=1,n_threads
#ifndef _CYLINDRICAL_
              uf1 => this%uf(:,lbound(this%uf,2)+this%m%jf(j-1):lbound(this%uf,2)+this%m%jf(j)-1)
              uf2 => source%uf(:,lbound(source%uf,2)+this%m%jf(j-1):lbound(source%uf,2)+this%m%jf(j)-1)
#else
              uf1 => this%uf(:,:,lbound(this%uf,3)+this%m%jf(j-1):lbound(this%uf,3)+this%m%jf(j)-1)
              uf2 => source%uf(:,:,lbound(source%uf,3)+source%m%jf(j-1):lbound(source%uf,3)+this%m%jf(j)-1)
#endif
              uf1 = uf2
             end do
!$OMP END PARALLEL DO 
#endif

#endif
        end if
        
        class default
           stop "E: wave functions not belonging to the same method"
        end select
    end subroutine copy

    
    function norm(this) result(n)
        class(_WF_), intent(inout) :: this
        real(kind=prec) :: n
        !!! TODO handle norm in frequency space without transforming
        call this%to_real_space 
#ifdef _REAL_        
        n = this%m%g%norm_real_gridfun(this%u)
#else
        n = this%m%g%norm_complex_gridfun(this%u)
#endif        
        
    end function norm


    function inner_product(this, wf) result(n)
        class(_WF_), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: wf 
        _COMPLEX_OR_REAL_(kind=prec) :: n

        select type (wf)
        class is (_WF_)
        if (.not.associated(wf%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if 
   
        call this%to_real_space 
        call wf%to_real_space 
#ifdef _REAL_        
        n = this%m%g%inner_product_real_gridfun(this%u, wf%u)
#else
        n = this%m%g%inner_product_complex_gridfun(this%u, wf%u)
#endif        
        class default
           stop "E: wave functions not belonging to the same method"
        end select
    end function inner_product



    function norm_in_frequency_space(this) result(n)
        class(_WF_), intent(inout) :: this
        real(kind=prec) :: n
        real(kind=prec), parameter :: pi = 3.1415926535897932384626433832795028841971693993751_prec
#ifdef _OPENMP
#ifndef _CYLINDRICAL_
#ifdef _REAL_
        complex(kind=prec), pointer :: ufc(:,:)
#else
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf(:,:)
#endif
#else
#ifdef _REAL_
        complex(kind=prec), pointer :: ufc(:,:,:)
#else
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf(:,:,:)
#endif
#endif
        integer :: j
#endif
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: n1 
#endif
        
        call this%to_frequency_space 

#ifdef _REAL_

         if (0>=this%m%nfthetamin) then !only relevant for MPI ...
#ifndef _CYLINDRICAL_
             this%ufc(:,0) = sqrt(0.5_prec)*this%ufc(:,0)
#else
             this%ufc(:,:,0) = sqrt(0.5_prec)*this%ufc(:,:,0)
#endif
         endif
#ifndef _OPENMP
#ifndef _CYLINDRICAL_
          n = sum(real(this%ufc(this%m%nfrmin:this%m%nfrmax,:),kind=prec)**2 &
                 +aimag(this%ufc(this%m%nfrmin:this%m%nfrmax,:))**2)
#else
          n = sum(real(this%ufc(this%m%nfrmin:this%m%nfrmax,:,:),kind=prec)**2 &
                 +aimag(this%ufc(this%m%nfrmin:this%m%nfrmax,:,:))**2)
#endif
#else
          n = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, ufc) REDUCTION(+:n)
          do j=1,n_threads
#ifndef _CYLINDRICAL_
              ufc => this%ufc(this%m%nfrmin:this%m%nfrmax,&
                          lbound(this%ufc,2)+this%m%jf(j-1):lbound(this%ufc,2)+this%m%jf(j)-1)
#else
              ufc => this%ufc(this%m%nfrmin:this%m%nfrmax,:,&
                          lbound(this%ufc,3)+this%m%jf(j-1):lbound(this%ufc,3)+this%m%jf(j)-1)
#endif
              n = n + sum(real(ufc,kind=prec)**2 +aimag(ufc)**2)
          end do
!$OMP END PARALLEL DO 
#endif
         if (0>=this%m%nfthetamin) then !only relevant for MPI ...
#ifndef _CYLINDRICAL_
             this%ufc(:,0) = sqrt(2.0_prec)*this%ufc(:,0)
#else
             this%ufc(:,:,0) = sqrt(2.0_prec)*this%ufc(:,:,0)
#endif
         endif
         n = 2.0_prec*n

#else

#ifndef _OPENMP
#ifndef _CYLINDRICAL_
             n = sum(real(this%uf(this%m%nfrmin:this%m%nfrmax,:),kind=prec)**2 &
                 +aimag(this%uf(this%m%nfrmin:this%m%nfrmax,:))**2)
#else
             n = sum(real(this%uf(this%m%nfrmin:this%m%nfrmax,:,:),kind=prec)**2 &
                 +aimag(this%uf(this%m%nfrmin:this%m%nfrmax,:,:))**2)
#endif
#else
             n = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, uf) REDUCTION(+:n) 
             do j=1,n_threads
#ifndef _CYLINDRICAL_
              uf => this%uf(this%m%nfrmin:this%m%nfrmax,&
                         lbound(this%uf,2)+this%m%jf(j-1):lbound(this%uf,2)+this%m%jf(j)-1)
#else
              uf => this%uf(this%m%nfrmin:this%m%nfrmax,:,&
                         lbound(this%uf,3)+this%m%jf(j-1):lbound(this%uf,3)+this%m%jf(j)-1)
#endif
              n = n + sum(real(uf,kind=prec)**2 +aimag(uf)**2)
             end do
!$OMP END PARALLEL DO 
#endif

#endif

#ifdef _MPI_
        call MPI_Reduce(n, n1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            n = sqrt(n1) 
        end if 
        call MPI_Bcast(N, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#else 
        n = sqrt(n)
#endif
        
    end function norm_in_frequency_space


    subroutine normalize(this, norm)
        class(_WF_), intent(inout) :: this
        real(kind=prec), intent(out), optional :: norm
        real(kind=prec) :: n

        n = this%norm()
        if (present(norm)) then
            norm = n
        end if    

#ifdef _REAL_        
        call this%scale(1.0_prec/n)
#else        
        call this%scale(cmplx(1.0_prec/n, 0.0_prec, kind=prec))
#endif            
    end subroutine normalize


    function distance(this, wf) result(n)
        class(_WF_), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: wf 
        real(kind=prec) :: n
        select type (wf)
        class is (_WF_)
        if (.not.associated(wf%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if  
        if (associated(this%u, wf%u)) then
            n = 0.0_prec
            return
        end if

#ifdef _REAL_
        call this%axpy(wf, -1.0_prec)
        n = this%norm()
        call this%axpy(wf, +1.0_prec)
#else
        call this%axpy(wf, cmplx(-1.0_prec, 0.0_prec, kind=prec))
        n = this%norm()
        call this%axpy(wf, cmplx(-1.0_prec, 0.0_prec, kind=prec))
#endif        

        class default
           stop "E: wave functions not belonging to the same method"
        end select
    end function distance


    subroutine axpy(this, other, factor)
         class(_WF_), intent(inout) :: this
         class(_WAVE_FUNCTION_), intent(inout) :: other
         _COMPLEX_OR_REAL_(kind=prec), intent(in) :: factor
#ifdef _OPENMP
#ifndef _CYLINDRICAL_
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u1(:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u2(:,:)
#else
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u1(:,:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u2(:,:,:)
#endif
        integer :: j
#endif        
        select type (other)
        class is (_WF_)
        if (.not.associated(other%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if 
        call this%to_real_space 
        call other%to_real_space 
#ifndef _OPENMP
#ifdef _REAL_
        if(factor==1.0_prec) then
           this%u = this%u + other%u
        elseif(factor==-1.0_prec) then
           this%u = this%u - other%u
        else
           this%u = this%u + factor*other%u
        end if   
#else
        if(aimag(factor)==0.0_prec) then
            if(real(factor,prec)==1.0_prec) then
               this%u = this%u + other%u
            elseif(real(factor,prec)==-1.0_prec) then
               this%u = this%u - other%u
            else
               this%u = this%u + real(factor, kind=prec)*other%u
            end if   
        else
            this%u = this%u + factor*other%u
        end if
#endif        
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
         do j=1,n_threads
#ifndef _CYLINDRICAL_
              u1 => this%u(:,lbound(this%u,2)+this%m%g%jj(j-1):lbound(this%u,2)+this%m%g%jj(j)-1)
              u2 => other%u(:,lbound(other%u,2)+other%m%g%jj(j-1):lbound(other%u,2)+other%m%g%jj(j)-1)
#else
              u1 => this%u(:,:,lbound(this%u,3)+this%m%g%jj(j-1):lbound(this%u,3)+this%m%g%jj(j)-1)
              u2 => other%u(:,:,lbound(other%u,3)+other%m%g%jj(j-1):lbound(other%u,3)+other%m%g%jj(j)-1)
#endif
#ifdef _REAL_
              if(factor==1.0_prec) then
                 u1 = u1 + u2
              elseif(factor==-1.0_prec) then
                 u1 = u1 - u2
              else
                 u1 = u1 + factor*u2
              end if   
#else
              if(aimag(factor)==0.0_prec) then
                  if(real(factor,prec)==1.0_prec) then
                     u1 = u1 + u2
                  elseif(real(factor,prec)==-1.0_prec) then
                     u1 = u1 - u2
                  else
                     u1 = u1 + real(factor, kind=prec)*u2
                  end if   
              else
                  u1 = u1 + factor*u2
              end if
#endif       
         end do
!$OMP END PARALLEL DO 
#endif
        class default
           stop "E: wave functions not belonging to the same method"
        end select
    end subroutine axpy
   

    subroutine scale(this, factor)
        class(_WF_), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: factor   
#ifdef _OPENMP
#ifndef _CYLINDRICAL_
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:)
#else        
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:,:)
#endif
        integer :: j
#endif        
#ifndef _OPENMP
#ifndef _REAL_
        if(aimag(factor)==0.0_prec) then
            this%u = real(factor, kind=prec)*this%u
        else
#endif        
            this%u = factor*this%u
#ifndef _REAL_
        end if
#endif        
#else
!$OMP PARALLEL DO PRIVATE(j, u) 
         do j=1,n_threads
#ifndef _CYLINDRICAL_
              u => this%u(:,lbound(this%u,2)+this%m%g%jj(j-1):lbound(this%u,2)+this%m%g%jj(j)-1)
#else
              u => this%u(:,:,lbound(this%u,3)+this%m%g%jj(j-1):lbound(this%u,3)+this%m%g%jj(j)-1)
#endif
#ifndef _REAL_
              if(aimag(factor)==0.0_prec) then
                  u = real(factor,kind=prec)*u
              else
#endif
                  u = factor*u
#ifndef _REAL_
              end if
#endif        
         end do
!$OMP END PARALLEL DO 
#endif
    end subroutine scale
    
end module _MODULE_ 
