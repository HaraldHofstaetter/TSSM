!
! cpp -D_DIM_=2 -D_REAL_ -P -std=c89  tssm_fourier_CPP.F90 
!
#ifdef _REAL_
#define S0(x,y)  x ## _real_ ## y ## D 
#else
#define S0(x,y)  x ## _ ## y ## D 
#endif
#define S1(x,y) S0(x,y)
#define S(x) S1(x,_DIM_)
#ifdef _REAL_
 #define _WAVE_FUNCTION_ real_wave_function
 #define _COMPLEX_OR_REAL_ real
#else
 #define _WAVE_FUNCTION_ wave_function
 #define _COMPLEX_OR_REAL_ complex
#endif

#ifdef _QUADPRECISION_
#define fftw_alloc_real fftwq_alloc_real
#define fftw_alloc_complex fftwq_alloc_complex
#define fftw_plan_dft_r2c_1d fftwq_plan_dft_r2c_1d
#define fftw_plan_dft_c2r_1d fftwq_plan_dft_c2r_1d
#define fftw_plan_dft_1d fftwq_plan_dft_1d
#define fftw_plan_dft_r2c_2d fftwq_plan_dft_r2c_2d
#define fftw_plan_dft_c2r_2d fftwq_plan_dft_c2r_2d
#define fftw_plan_dft_2d fftwq_plan_dft_2d
#define fftw_plan_dft_r2c_3d fftwq_plan_dft_r2c_3d
#define fftw_plan_dft_c2r_3d fftwq_plan_dft_c2r_3d
#define fftw_plan_dft_3d fftwq_plan_dft_3d
#define fftw_plan_many_r2r fftwq_plan_many_r2r
#define fftw_execute_dft_r2c fftwq_execute_dft_r2c
#define fftw_execute_dft_c2r fftwq_execute_dft_c2r
#define fftw_execute_dft fftwq_execute_dft
#define fftw_execute_r2r fftwq_execute_r2r
#define fftw_destroy_plan fftwq_destroy_plan
#define fftw_free fftwq_free
#endif 

#ifdef _QUADPRECISION_
module S(tssmq_fourier)
    use tssmq_base, only: spectral_method, _WAVE_FUNCTION_
    use tssmq_grid
    use tssmq_fourier_common
#else
module S(tssm_fourier)
    use tssm_base, only: spectral_method, _WAVE_FUNCTION_
    use tssm_grid
    use tssm_fourier_common
#endif    
    use, intrinsic :: iso_c_binding, only: c_ptr
    implicit none

    private
    public :: S(fourier), S(wf_fourier)


    type, extends(spectral_method) :: S(fourier)
#if(_DIM_==1)
        type(grid_equidistant_1D) :: g 
#elif(_DIM_==2)
        type(grid_equidistant_2D) :: g 
#elif(_DIM_==3)
        type(grid_equidistant_3D) :: g 
#endif 
        integer :: nf1min !integer range for transformed data
        integer :: nf1max
#if(_DIM_>=2)
        integer :: nf2min
        integer :: nf2max
#endif
#if(_DIM_>=3)
        integer :: nf3min
        integer :: nf3max
#endif
!        integer(kind=C_SIZE_T) :: alloc_size

        integer :: boundary_conditions = periodic ! default dft

        real(kind=prec), pointer :: eigenvalues1(:) => null()
#if(_DIM_>=2)
        real(kind=prec), pointer :: eigenvalues2(:) => null()
#endif 
#if(_DIM_>=3)
        real(kind=prec), pointer :: eigenvalues3(:) => null()
#endif 

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
        procedure :: finalize => finalize_method
        !final :: final_fourier_1D
        !! Fortran 2003 feature final seems to be not properly implemented
        !! in the gcc/gfortran compiler :(
    end type S(fourier)

    interface  S(fourier) ! constructor
        module procedure new_method
    end interface S(fourier)


    type, extends(_WAVE_FUNCTION_) :: S(wf_fourier)
        class(S(fourier)), pointer :: m 

#ifdef _REAL_
        real(kind=prec), pointer  :: up(:)
        complex(kind=prec), pointer  :: ucp(:)
#else
        complex(kind=prec), pointer  :: up(:)
        real(kind=prec), pointer  :: urp(:)
#endif

#if(_DIM_==1)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf(:)
#elif(_DIM_==2)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf(:,:)
#elif(_DIM_==3)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf(:,:,:)
#endif

#ifdef _REAL_
#if(_DIM_==1)
        complex(kind=prec), pointer :: ufc(:)
#elif(_DIM_==2)
        complex(kind=prec), pointer :: ufc(:,:)
#elif(_DIM_==3)
        complex(kind=prec), pointer :: ufc(:,:,:)
#endif
#endif

        type(c_ptr) :: plan_forward 
        type(c_ptr) :: plan_backward

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
        procedure :: distance
        procedure :: axpy
        procedure :: scale
    end type S(wf_fourier)

    interface S(wf_fourier) ! constructor
        module procedure new_wf
    end interface S(wf_fourier)


contains

#if(_DIM_==1)
    function new_method(nx, xmin, xmax, boundary_conditions) result(this)
#elif(_DIM_==2)
    function new_method(nx, xmin, xmax, ny, ymin, ymax, boundary_conditions) result(this)
#elif(_DIM_==3)
    function new_method(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax, boundary_conditions) result(this)
#endif
        type(S(fourier)) :: this
        real(kind=prec),  intent(in) :: xmin 
        real(kind=prec),  intent(in) :: xmax
        integer, intent(in) :: nx
#if(_DIM_>=2)
        real(kind=prec),  intent(in) :: ymin 
        real(kind=prec),  intent(in) :: ymax
        integer, intent(in) :: ny
#endif
#if(_DIM_>=3)
        real(kind=prec),  intent(in) :: zmin 
        real(kind=prec),  intent(in) :: zmax
        integer, intent(in) :: nz
#endif
        integer, intent(in), optional :: boundary_conditions
        integer :: i

#ifdef _MPI_
#if(_DIM_==1)
        integer(C_SIZE_T) :: local_ni, local_i_start
        integer(C_SIZE_T) :: local_no, local_o_start
#else
        integer(C_SIZE_T) :: local_n_last, local_last_start
        integer(C_SIZE_T) :: local_n_2ndlast_trans, local_2ndlast_start_trans
#endif
#endif
        if (present(boundary_conditions)) then
            this%boundary_conditions = boundary_conditions
        end if
        
#if(_DIM_==1)
        this%g = grid_equidistant_1D(nx, xmin, xmax)
#elif(_DIM_==2)
        this%g = grid_equidistant_2D(nx, xmin, xmax, ny, ymin, ymax)
#elif(_DIM_==3)
        this%g = grid_equidistant_3D(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax)
#endif

        select case(this%boundary_conditions) ! adjust grid
        case (periodic)
            this%g%nn1min = 1
            this%g%nn1max = this%g%nx
#if(_DIM_>=2)
            this%g%nn2min = 1
            this%g%nn2max = this%g%ny
#endif
#if(_DIM_>=3)
            this%g%nn3min = 1
            this%g%nn3max = this%g%nz
#endif            
        case (dirichlet)
            this%g%nn1min = 1
            this%g%nn1max = this%g%nx-1
#if(_DIM_>=2)
            this%g%nn2min = 1
            this%g%nn2max = this%g%ny-1
#endif
#if(_DIM_>=3)
            this%g%nn3min = 1
            this%g%nn3max = this%g%nz-1
#endif
        case (neumann)
            this%g%nn1min = 0
            this%g%nn1max = this%g%nx
#if(_DIM_>=2)
            this%g%nn2min = 0
            this%g%nn2max = this%g%ny
#endif
#if(_DIM_>=3)
            this%g%nn3min = 0
            this%g%nn3max = this%g%nz
#endif
        end select


#ifdef _MPI_
#ifdef _REAL_
        select case(this%boundary_conditions) ! adjust grid
        case (periodic)
#if(_DIM_==1)

#elif(_DIM_==2)
        this%g%n1min = this%g%nn1min
        this%g%n1max = this%g%nn1max
        this%g%alloc_size = 2*fftw_mpi_local_size_2d_transposed( &
               int(this%g%ny, kind=C_SIZE_T), &
               int(this%g%nx/2+1, kind=C_SIZE_T), &
               MPI_COMM_WORLD, &
               local_n_last, local_last_start, &
               local_n_2ndlast_trans, local_2ndlast_start_trans)
        this%nf1min = this%g%nn2min ! 2->1 ... transposed
        this%nf1max = this%g%nn2max 
        this%g%n2min = local_last_start + this%g%n2min
        this%g%n2max = this%g%n2min + local_n_last - 1
        this%nf2min = local_2ndlast_start_trans + this%g%n1min
        this%nf2max = this%nf2min + local_n_2ndlast_trans - 1

        this%g%m1min = this%g%n1min
        this%g%m1max = this%g%n1max+2
        this%g%m2min = this%g%n2min
        this%g%m2max = this%g%n2max
#elif(_DIM_==3)
        this%g%n1min = this%g%nn1min
        this%g%n1max = this%g%nn1max
        this%g%n2min = this%g%nn2min
        this%g%n2max = this%g%nn2max
        this%g%n3min = this%g%nn3min
        this%g%n3max = this%g%nn3max
        this%g%alloc_size = 2*fftw_mpi_local_size_3d_transposed( &
               int(this%g%nz, kind=C_SIZE_T), &
               int(this%g%ny, kind=C_SIZE_T), &
               int(this%g%nx/2+1, kind=C_SIZE_T), &
               MPI_COMM_WORLD, &
               local_n_last, local_last_start, &
               local_n_2ndlast_trans, local_2ndlast_start_trans)
        this%nf1min = this%g%nn1min 
        this%nf1max = this%g%nx/2+1
        this%nf2min = this%g%nn3min ! 3->2 ... transposed
        this%nf2max = this%g%nn3max 
        this%g%n3min = local_last_start + this%g%n3min
        this%g%n3max = this%g%n3min + local_n_last - 1
        this%nf3min = local_2ndlast_start_trans + this%g%n2min
        this%nf3max = this%nf3min + local_n_2ndlast_trans - 1

        this%g%m1min = this%g%n1min
        this%g%m1max = this%g%n1max+2
        this%g%m2min = this%g%n2min
        this%g%m2max = this%g%n2max
        this%g%m3min = this%g%n3min
        this%g%m3max = this%g%n3max
#endif 



        case (dirichlet, neumann)
#endif 
#if(_DIM_==1)
!TODO: Check that dimensions and n_proc are compatible !!!!
        this%g%alloc_size = fftw_mpi_local_size_1d( &
               int(this%g%nn1max-this%g%nn1min+1, kind=C_SIZE_T), &
               MPI_COMM_WORLD, &
               FFTW_FORWARD, FFTW_MPI_SCRAMBLED_IN, &
               local_ni, local_i_start, &
               local_no ,local_o_start ) 
        this%g%n1min = local_i_start + 1
        this%g%n1max = local_i_start + local_ni 
        this%nf1min = local_o_start + 1
        this%nf1max = local_o_start + local_no
        this%g%m1min = this%g%n1min
        this%g%m1max = this%g%n1max
#elif(_DIM_==2)
        this%g%n1min = this%g%nn1min
        this%g%n1max = this%g%nn1max
        this%g%n2min = this%g%nn2min
        this%g%n2max = this%g%nn2max
        this%g%alloc_size = fftw_mpi_local_size_2d_transposed( &
               int(this%g%nn2max-this%g%nn2min+1, kind=C_SIZE_T), &
               int(this%g%nn1max-this%g%nn1min+1, kind=C_SIZE_T), &
               MPI_COMM_WORLD, &
               local_n_last, local_last_start, &
               local_n_2ndlast_trans, local_2ndlast_start_trans)
        this%nf1min = this%g%nn2min ! 2->1 ... transposed
        this%nf1max = this%g%nn2max 
        this%g%n2min = local_last_start + this%g%n2min
        this%g%n2max = this%g%n2min + local_n_last - 1
        this%nf2min = local_2ndlast_start_trans + this%g%n1min
        this%nf2max = this%nf2min + local_n_2ndlast_trans - 1

        this%g%m2min = this%g%n2min
        this%g%m2max = this%g%n2max
!print *, "P", this_proc, this%g%alloc_size, this%g%n1min, this%g%n1max, this%g%n2min, this%g%n2max, &
!              this%nf1min, this%nf1max, this%nf2min, this%nf2max 
#elif(_DIM_==3)
        this%g%n1min = this%g%nn1min
        this%g%n1max = this%g%nn1max
        this%g%n2min = this%g%nn2min
        this%g%n2max = this%g%nn2max
        this%g%n3min = this%g%nn3min
        this%g%n3max = this%g%nn3max
        this%g%alloc_size = fftw_mpi_local_size_3d_transposed( &
               int(this%g%nn3max-this%g%nn3min+1, kind=C_SIZE_T), &
               int(this%g%nn2max-this%g%nn2min+1, kind=C_SIZE_T), &
               int(this%g%nn1max-this%g%nn1min+1, kind=C_SIZE_T), &
               MPI_COMM_WORLD, &
               local_n_last, local_last_start, &
               local_n_2ndlast_trans, local_2ndlast_start_trans)
        this%nf1min = this%g%nn1min 
        this%nf1max = this%g%nn1max 
        this%nf2min = this%g%nn3min ! 3->2 ... transposed
        this%nf2max = this%g%nn3max 
        this%g%n3min = local_last_start + this%g%n3min
        this%g%n3max = this%g%n3min + local_n_last - 1
        this%nf3min = local_2ndlast_start_trans + this%g%n2min
        this%nf3max = this%nf3min + local_n_2ndlast_trans - 1

        this%g%m1min = this%g%n1min
        this%g%m1max = this%g%n1max
        this%g%m2min = this%g%n2min
        this%g%m2max = this%g%n2max
        this%g%m3min = this%g%n3min
        this%g%m3max = this%g%n3max
#endif 
#ifdef _REAL_
        end select
#endif 

#else
        this%g%n1min = this%g%nn1min
        this%g%n1max = this%g%nn1max
        this%g%m1min = this%g%nn1min
        this%nf1min = this%g%n1min
#ifdef _REAL_
        select case(this%boundary_conditions) ! adjust grid
        case (periodic)
            this%g%m1max =  2*(this%g%nx/2+1)  ! for indexing REAL data
            this%nf1max = (this%g%nx/2+1)      ! for indexing complex data
        case (dirichlet, neumann)
#endif 
            this%g%m1max = this%g%nn1max
            this%nf1max = this%g%n1max
#ifdef _REAL_
        end select
#endif 

#if(_DIM_>=2)
        this%g%n2min = this%g%nn2min
        this%g%n2max = this%g%nn2max
        this%g%m2min = this%g%nn2min
        this%g%m2max = this%g%nn2max
        this%nf2min = this%g%n2min
        this%nf2max = this%g%n2max
#endif
#if(_DIM_>=3)
        this%g%n3min = this%g%nn3min
        this%g%n3max = this%g%nn3max
        this%g%m3min = this%g%nn3min
        this%g%m3max = this%g%nn3max
        this%nf3min = this%g%n3min
        this%nf3max = this%g%n3max
#endif

#if(_DIM_==1)
        this%g%alloc_size = this%g%m1max-this%g%m1min+1
#elif(_DIM_==2)
        this%g%alloc_size = (this%g%m1max-this%g%m1min+1)*(this%g%m2max-this%g%m2min+1)
#elif(_DIM_==3)
        this%g%alloc_size = (this%g%m1max-this%g%m1min+1)*(this%g%m2max-this%g%m2min+1) &
                         *(this%g%m3max-this%g%m3min+1)
#endif
#endif

        allocate( this%g%nodes_x(this%g%n1min:this%g%n1max ) )
        this%g%nodes_x = (/ ( this%g%xmin+this%g%dx*real(i, prec), i = this%g%n1min, this%g%n1max )  /)
#if(_DIM_>=2)
        allocate( this%g%nodes_y(this%g%n2min:this%g%n2max ) )
        this%g%nodes_y = (/ ( this%g%ymin+this%g%dy*real(i, prec), i = this%g%n2min, this%g%n2max )  /)
#endif
#if(_DIM_>=3)
        allocate( this%g%nodes_z(this%g%n3min:this%g%n3max ) )
        this%g%nodes_z = (/ ( this%g%zmin+this%g%dz*real(i, prec), i = this%g%n3min, this%g%n3max )  /)
#endif

        allocate( this%eigenvalues1(this%nf1min:this%nf1max ) ) 
#if (defined(_MPI_)&&(_DIM_==1))       
        call get_scrambled_eigenvalues(this%eigenvalues1, this%g%xmax-this%g%xmin, &
                             this%g%nx, this%nf1min, this%nf1max, &
                             this%g%alloc_size)
#elif (defined(_MPI_)&&(_DIM_==2))       
        call get_eigenvalues(this%eigenvalues1, this%g%ymax-this%g%ymin, &
                             this%g%ny, this%nf1min, this%nf1max, &   ! 2<->1 transposition
                             this%boundary_conditions)
#else                             
        call get_eigenvalues(this%eigenvalues1, this%g%xmax-this%g%xmin, &
                             this%g%nx, this%nf1min, this%nf1max, &
                             this%boundary_conditions)
#endif                              
#if(_DIM_>=2)
        allocate( this%eigenvalues2(this%nf2min:this%nf2max ) ) 
#if (defined(_MPI_)&&(_DIM_==2))        
        call get_eigenvalues(this%eigenvalues2, this%g%xmax-this%g%xmin, & ! 2<->1 transposition
                             this%g%nx, this%nf2min, this%nf2max, &
                             this%boundary_conditions)
#elif (defined(_MPI_)&&(_DIM_==3))        
        call get_eigenvalues(this%eigenvalues2, this%g%zmax-this%g%zmin, & ! 3<->2 transposition
                             this%g%nz, this%nf2min, this%nf2max, &
                             this%boundary_conditions)
#else         
        call get_eigenvalues(this%eigenvalues2, this%g%ymax-this%g%ymin, &
                             this%g%ny, this%nf2min, this%nf2max, &
                             this%boundary_conditions)
#endif                              
#endif
#if(_DIM_>=3)
        allocate( this%eigenvalues3(this%nf3min:this%nf3max ) ) 
#ifdef _MPI_
        call get_eigenvalues(this%eigenvalues3, this%g%ymax-this%g%ymin, & ! 3<->2 transposition
                             this%g%ny, this%nf3min, this%nf3max, &
                             this%boundary_conditions)
#else
        call get_eigenvalues(this%eigenvalues3, this%g%zmax-this%g%zmin, &
                             this%g%nz, this%nf3min, this%nf3max, &
                             this%boundary_conditions)
#endif
#endif
#ifdef _OPENMP
        allocate( this%g%jj(0:n_threads) )
        allocate( this%jf(0:n_threads) )
#if(_DIM_==1)
        do i=0,n_threads-1
           this%g%jj(i) = i*ceiling(real(this%g%n1max-this%g%n1min+1)/real(n_threads))
           this%jf(i) = i*ceiling(real(this%nf1max-this%nf1min+1)/real(n_threads))
        end do
        this%g%jj(n_threads) = this%g%n1max-this%g%n1min+1 
        this%jf(n_threads) = this%nf1max-this%nf1min+1 
#elif(_DIM_==2)
        do i=0,n_threads-1
           this%g%jj(i) = i*ceiling(real(this%g%n2max-this%g%n2min+1)/real(n_threads))
           this%jf(i) = i*ceiling(real(this%nf2max-this%nf2min+1)/real(n_threads))
        end do
        this%g%jj(n_threads) = this%g%n2max-this%g%n2min+1
        this%jf(n_threads) = this%nf2max-this%nf2min+1 
#elif(_DIM_==3)
        do i=0,n_threads-1
           this%g%jj(i) = i*ceiling(real(this%g%n3max-this%g%n3min+1)/real(n_threads))
           this%jf(i) = i*ceiling(real(this%nf3max-this%nf3min+1)/real(n_threads))
        end do
        this%g%jj(n_threads) = this%g%n3max-this%g%n3min+1
        this%jf(n_threads) = this%nf3max-this%nf3min+1 
#endif

#endif
    end function new_method
 

    subroutine finalize_method(this)
        class(S(fourier)), intent(inout) :: this

        deallocate( this%g%nodes_x )
        deallocate( this%eigenvalues1 )
#if(_DIM_>=2)
        deallocate( this%g%nodes_y )
        deallocate( this%eigenvalues2 )
#endif
#if(_DIM_>=3)
        deallocate( this%g%nodes_z )
        deallocate( this%eigenvalues3 )
#endif
#ifdef _OPENMP
        deallocate( this%g%jj )
        deallocate( this%jf )
#endif        
    end subroutine finalize_method

    function new_wf(m, u, coefficient) result(this)
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_ptr, c_loc
        type(S(wf_fourier)) :: this
        class(S(fourier)), target, intent(inout) :: m
#ifdef _REAL_       
#if(_DIM_==1)
        real(kind=prec), optional, target, intent(inout) :: u(:)
#elif(_DIM_==2)
        real(kind=prec), optional, target, intent(inout) :: u(:,:)
#elif(_DIM_==3)
        real(kind=prec), optional, target, intent(inout) :: u(:,:,:)
#endif
        real(kind=prec), optional, intent(in) :: coefficient
#else            
#if(_DIM_==1)
        complex(kind=prec), optional, target, intent(inout) :: u(:)
#elif(_DIM_==2)
        complex(kind=prec), optional, target, intent(inout) :: u(:,:)
#elif(_DIM_==3)
        complex(kind=prec), optional, target, intent(inout) :: u(:,:,:)
#endif
        complex(kind=prec), optional, intent(in) :: coefficient
#endif
#ifdef _REAL_       
#if(_DIM_==1)
        real(kind=prec), pointer :: umem(:)
#elif(_DIM_==2)
        real(kind=prec), pointer  :: umem(:,:)
#elif(_DIM_==3)
        real(kind=prec), pointer  :: umem(:,:,:)
#endif
#else            
#if(_DIM_==1)
        complex(kind=prec), pointer  :: umem(:)
#elif(_DIM_==2)
        complex(kind=prec), pointer  :: umem(:,:)
#elif(_DIM_==3)
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
#if(_DIM_==1)
            umem => u(m%g%m1min:m%g%m1max) 
#elif(_DIM_==2)
            umem => u(m%g%m1min:m%g%m1max, m%g%m2min:m%g%m2max) 
#elif(_DIM_==3)
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
#if(_DIM_==1)
            umem(this%m%g%m1min:this%m%g%m1max) => this%up 
            this%u => umem(this%m%g%n1min:this%m%g%n1max) 
#elif(_DIM_==2)
            umem(this%m%g%m1min:this%m%g%m1max, this%m%g%m2min:this%m%g%m2max) => this%up
            this%u => umem(this%m%g%n1min:this%m%g%n1max, this%m%g%n2min:this%m%g%n2max)
#elif(_DIM_==3)
            umem(this%m%g%m1min:this%m%g%m1max, this%m%g%m2min:this%m%g%m2max, this%m%g%m3min:this%m%g%m3max) => this%up 
            this%u => umem(this%m%g%n1min:this%m%g%n1max, this%m%g%n2min:this%m%g%n2max, this%m%g%n3min:this%m%g%n3max)  
#endif
        end if

#if(_DIM_==1)
        this%uf(this%m%nf1min:this%m%nf1max) => this%up 
#elif(_DIM_==2)
        this%uf(this%m%nf1min:this%m%nf1max, this%m%nf2min:this%m%nf2max) => this%up
#elif(_DIM_==3)
        this%uf(this%m%nf1min:this%m%nf1max, this%m%nf2min:this%m%nf2max, this%m%nf3min:this%m%nf3max) =>this%up  
#endif
#ifdef _REAL_
        select case(this%m%boundary_conditions) 
        case (periodic)
#if(_DIM_==1)
        this%ufc(this%m%nf1min:this%m%nf1max) => this%ucp 
#elif(_DIM_==2)
        this%ufc(this%m%nf1min:this%m%nf1max, this%m%nf2min:this%m%nf2max) => this%ucp 
#elif(_DIM_==3)
        this%ufc(this%m%nf1min:this%m%nf1max, this%m%nf2min:this%m%nf2max, this%m%nf3min:this%m%nf3max) => this%ucp 
#endif 
        case (dirichlet, neumann)
        case default
        end select
#endif 
        call this%create_plans
    end function new_wf


    function clone(this) 
        class(S(wf_fourier)), intent(inout) :: this
        class(_WAVE_FUNCTION_), pointer :: clone
        type(S(wf_fourier)), pointer :: p

        allocate( p )
        p = S(wf_fourier)(this%m, coefficient=this%coefficient)
        clone => p
    end function clone


    subroutine finalize_wf(this)
        use, intrinsic :: iso_c_binding, only: c_loc, c_associated
        class(S(wf_fourier)), intent(inout) :: this

        if (.not.c_associated(this%plan_backward, this%plan_forward)) then 
            call fftw_destroy_plan(this%plan_backward)
        end if
        call fftw_destroy_plan(this%plan_forward)
        !call fftw_free(c_loc(this%up))
        call fftw_free(c_loc(this%up(1)))
    end subroutine finalize_wf

    
    subroutine create_plans(this)
        class(S(wf_fourier)) :: this
#ifdef _MPI_
        integer(kind=C_SIZE_T) :: dft_dim(_DIM_)  
#else
        integer :: dft_dim(_DIM_)  
#endif
        integer :: dft_kind(_DIM_) 

        call load_fftw_wisdom

        select case(this%m%boundary_conditions) 
        case (periodic)
#if(_DIM_==1)
#ifdef _REAL_
            this%plan_forward = fftw_plan_dft_r2c_1d(this%m%g%nx, this%up, this%ucp, fftw_planning_rigor)
            this%plan_backward = fftw_plan_dft_c2r_1d(this%m%g%nx, this%ucp, this%up, fftw_planning_rigor)

#else
#ifdef _MPI_
            this%plan_forward = fftw_mpi_plan_dft_1d(int(this%m%g%nx, kind=C_SIZE_T), this%up, this%up, &
                             MPI_COMM_WORLD, FFTW_FORWARD, fftw_planning_rigor+FFTW_MPI_SCRAMBLED_OUT)
            this%plan_backward = fftw_mpi_plan_dft_1d(int(this%m%g%nx, kind=C_SIZE_T), this%up, this%up, &
                             MPI_COMM_WORLD, FFTW_BACKWARD, fftw_planning_rigor+FFTW_MPI_SCRAMBLED_IN)
#else
            this%plan_forward = fftw_plan_dft_1d(this%m%g%nx, this%up, this%up, FFTW_FORWARD, fftw_planning_rigor)
            this%plan_backward = fftw_plan_dft_1d(this%m%g%nx, this%up, this%up, FFTW_BACKWARD, fftw_planning_rigor)
#endif            
#endif            
#elif(_DIM_==2)
#ifdef _REAL_
#ifdef _MPI_
            this%plan_forward = fftw_mpi_plan_dft_r2c_2d(int(this%m%g%ny, kind=C_SIZE_T), &
                                      int(this%m%g%nx, kind=C_SIZE_T), this%up, this%ucp, &
                                      MPI_COMM_WORLD, fftw_planning_rigor+FFTW_MPI_TRANSPOSED_OUT)
            this%plan_backward = fftw_mpi_plan_dft_c2r_2d(int(this%m%g%ny, kind=C_SIZE_T), &
                                      int(this%m%g%nx, kind=C_SIZE_T), this%ucp, this%up,  &
                                      MPI_COMM_WORLD, fftw_planning_rigor+FFTW_MPI_TRANSPOSED_IN)
#else
            this%plan_forward = fftw_plan_dft_r2c_2d(this%m%g%ny, this%m%g%nx, this%up, this%ucp, &
                                                     fftw_planning_rigor)
            this%plan_backward = fftw_plan_dft_c2r_2d(this%m%g%ny, this%m%g%nx, this%ucp, this%up,&
                                                      fftw_planning_rigor)
#endif                                                      
#else
            !Note the reversed order ny, nx. For this cf. FFTW docu/section 7.2
#ifdef _MPI_
            this%plan_forward = fftw_mpi_plan_dft_2d(int(this%m%g%ny, kind=C_SIZE_T), &
                 int(this%m%g%nx, kind=C_SIZE_T), this%up, this%up, MPI_COMM_WORLD, &
                 FFTW_FORWARD, fftw_planning_rigor + FFTW_MPI_TRANSPOSED_OUT)
            this%plan_backward = fftw_mpi_plan_dft_2d(int(this%m%g%ny, kind=C_SIZE_T), &
                 int(this%m%g%nx, kind=C_SIZE_T), this%up, this%up, MPI_COMM_WORLD, &
                 FFTW_BACKWARD, fftw_planning_rigor + FFTW_MPI_TRANSPOSED_IN)
#else
            this%plan_forward = fftw_plan_dft_2d(this%m%g%ny, this%m%g%nx, &
                                                 this%up, this%up, FFTW_FORWARD, fftw_planning_rigor)
            this%plan_backward = fftw_plan_dft_2d(this%m%g%ny, this%m%g%nx,&
                                                  this%up, this%up, FFTW_BACKWARD, fftw_planning_rigor)
#endif                                                  
#endif            
#elif(_DIM_==3)
#ifdef _REAL_
#ifdef _MPI_
            this%plan_forward = fftw_mpi_plan_dft_r2c_3d(int(this%m%g%nz, kind=C_SIZE_T),   &
                                       int(this%m%g%ny, kind=C_SIZE_T), &
                                       int(this%m%g%nx, kind=C_SIZE_T), this%up, this%ucp, &
                                       MPI_COMM_WORLD, fftw_planning_rigor+FFTW_MPI_TRANSPOSED_OUT)
            this%plan_backward = fftw_mpi_plan_dft_c2r_3d(int(this%m%g%nz, kind=C_SIZE_T), &
                                       int(this%m%g%ny, kind=C_SIZE_T), &
                                       int(this%m%g%nx, kind=C_SIZE_T), this%ucp, this%up,  &
                                       MPI_COMM_WORLD, fftw_planning_rigor+FFTW_MPI_TRANSPOSED_IN)

#else
            this%plan_forward = fftw_plan_dft_r2c_3d(this%m%g%nz, this%m%g%ny, this%m%g%nx, this%up, &
                                         this%ucp,  fftw_planning_rigor)
            this%plan_backward = fftw_plan_dft_c2r_3d(this%m%g%nz, this%m%g%ny, this%m%g%nx, this%ucp,&
                                         this%up, fftw_planning_rigor)
#endif                                                      
#else
            !Note the reversed order nz, ny, nx. For this cf. FFTW docu/section 7.2
#ifdef _MPI_
            this%plan_forward = fftw_mpi_plan_dft_3d(int(this%m%g%nz, kind=C_SIZE_T), &
                 int(this%m%g%ny, kind=C_SIZE_T), int(this%m%g%nx, kind=C_SIZE_T),this%up, this%up, MPI_COMM_WORLD, &
                 FFTW_FORWARD, fftw_planning_rigor + FFTW_MPI_TRANSPOSED_OUT)
            this%plan_backward = fftw_mpi_plan_dft_3d(int(this%m%g%nz, kind=C_SIZE_T), &
                 int(this%m%g%ny, kind=C_SIZE_T), int(this%m%g%nx, kind=C_SIZE_T),this%up, this%up, MPI_COMM_WORLD, &
                 FFTW_BACKWARD, fftw_planning_rigor + FFTW_MPI_TRANSPOSED_IN)

#else
            this%plan_forward = fftw_plan_dft_3d(this%m%g%nz, this%m%g%ny, this%m%g%nx, &
                                                 this%up, this%up, FFTW_FORWARD, fftw_planning_rigor)
            this%plan_backward = fftw_plan_dft_3d(this%m%g%nz, this%m%g%ny, this%m%g%nx,&
                                                  this%up, this%up, FFTW_BACKWARD, fftw_planning_rigor)
#endif   
#endif   
#endif
        case (dirichlet, neumann)
#if(_DIM_==1)
            select case(this%m%boundary_conditions) 
                case (dirichlet)
                     dft_dim  =  (/ this%m%g%nx-1 /)
                     dft_kind = (/ FFTW_RODFT00 /)   !(O)dd
                case (neumann)
                     dft_dim  =  (/ this%m%g%nx+1 /)
                     dft_kind = (/ FFTW_REDFT00  /)   !(E)ven
            end select
#elif(_DIM_==2)
            select case(this%m%boundary_conditions) 
                case (dirichlet)
                     dft_dim  =  (/ this%m%g%ny-1, this%m%g%nx-1 /)
                     dft_kind = (/ FFTW_RODFT00, FFTW_RODFT00 /)   !(O)dd
                case (neumann)
                     dft_dim  =  (/ this%m%g%ny+1, this%m%g%nx+1 /)
                     dft_kind = (/ FFTW_REDFT00, FFTW_REDFT00 /)   !(E)ven
            end select
#elif(_DIM_==3)
            select case(this%m%boundary_conditions) 
                case (dirichlet)
                     dft_dim  =  (/ this%m%g%nz-1, this%m%g%ny-1, this%m%g%nx-1 /)
                     dft_kind = (/ FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00 /)   !(O)dd
                case (neumann)
                     dft_dim  =  (/ this%m%g%nz+1, this%m%g%ny+1, this%m%g%nx+1 /)
                     dft_kind = (/ FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00 /)   !(E)ven
            end select
#endif
#ifdef _REAL_
#ifdef _MPI_
            this%plan_forward = fftw_mpi_plan_many_r2r( _DIM_, dft_dim, 1_C_SIZE_T, &
                                     FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &
                                     this%up, this%up, MPI_COMM_WORLD, dft_kind, &
                                     fftw_planning_rigor + FFTW_MPI_TRANSPOSED_OUT)
            this%plan_backward = fftw_mpi_plan_many_r2r( _DIM_, dft_dim, 1_C_SIZE_T, &
                                     FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &
                                     this%up, this%up, MPI_COMM_WORLD, dft_kind, &
                                     fftw_planning_rigor + FFTW_MPI_TRANSPOSED_IN)

#else
            this%plan_forward = fftw_plan_many_r2r(_DIM_, dft_dim, 1, &
                                     this%up, dft_dim, 1, 1, &
                                     this%up, dft_dim, 1, 1, &
                                     dft_kind , fftw_planning_rigor)
            this%plan_backward = this%plan_forward 
#endif                                     
#else
#ifdef _MPI_
            this%plan_forward = fftw_mpi_plan_many_r2r( _DIM_, dft_dim, 2_C_SIZE_T, &
                                     FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &
                                     this%urp, this%urp, MPI_COMM_WORLD, dft_kind, &
                                     fftw_planning_rigor + FFTW_MPI_TRANSPOSED_OUT)
            this%plan_backward = fftw_mpi_plan_many_r2r( _DIM_, dft_dim, 2_C_SIZE_T, &
                                     FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &
                                     this%urp, this%urp, MPI_COMM_WORLD, dft_kind, &
                                     fftw_planning_rigor + FFTW_MPI_TRANSPOSED_IN)
#else
            this%plan_forward = fftw_plan_many_r2r(_DIM_, dft_dim, 2, &
                                     this%urp, dft_dim, 2, 1, &
                                     this%urp, dft_dim, 2, 1, &
                                     dft_kind , fftw_planning_rigor)
            this%plan_backward = this%plan_forward 
#endif
#endif
        end select

        call save_fftw_wisdom

    end subroutine create_plans


    subroutine to_frequency_space(this)
        class(S(wf_fourier)), intent(inout) :: this
        
        real(kind=prec), parameter :: f = sqrt(0.5_prec)

        if (this%is_real_space) then
            select case(this%m%boundary_conditions) 
            case (periodic)
#ifdef _MPI_
#ifdef _REAL_
                call fftw_mpi_execute_dft_r2c(this%plan_forward, this%up, this%ucp)
#else
                call fftw_mpi_execute_dft(this%plan_forward, this%up, this%up)
#endif
            case (dirichlet, neumann)
#ifdef _REAL_
                call fftw_mpi_execute_r2r(this%plan_forward, this%up, this%up)
#else
                call fftw_mpi_execute_r2r(this%plan_forward, this%urp, this%urp)
#endif
#else
#ifdef _REAL_
                call fftw_execute_dft_r2c(this%plan_forward, this%up, this%ucp)
#else
                call fftw_execute_dft(this%plan_forward, this%up, this%up)
#endif
            case (dirichlet, neumann)
#ifdef _REAL_
                call fftw_execute_r2r(this%plan_forward, this%up, this%up)
#else
                call fftw_execute_r2r(this%plan_forward, this%urp, this%urp)
#endif
#endif
            end select

            select case(this%m%boundary_conditions) 
            case (neumann)
#if(_DIM_==1)            
                if (this%m%nf1min==0) then
                     this%uf(0) = f * this%uf(0)
               end if      
#elif(_DIM_==2)          
               if (this%m%nf1min==0) then
                     this%uf(0,:) = f * this%uf(0,:)
               end if      
               if (this%m%nf2min==0) then
                     this%uf(:,0) = f * this%uf(:,0)
               end if      
#elif(_DIM_==3)
               if (this%m%nf1min==0) then
                     this%uf(0,:,:) = f * this%uf(0,:,:)
               end if      
               if (this%m%nf2min==0) then
                     this%uf(:,0,:) = f * this%uf(:,0,:)
               end if      
               if (this%m%nf3min==0) then
                     this%uf(:,:,0) = f * this%uf(:,:,0)
               end if     
#endif
            end select

            this%is_real_space = .false.
        end if    
    end subroutine to_frequency_space



    subroutine to_real_space(this)
        class(S(wf_fourier)), intent(inout) :: this

        real(kind=prec), parameter :: f = sqrt(2.0_prec)

#ifdef _OPENMP
#if(_DIM_==1)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:)
#elif(_DIM_==2)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:)
#elif(_DIM_==3)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:,:)
#endif
        integer :: j
#endif        
        
        
        if (.not. this%is_real_space) then

            select case(this%m%boundary_conditions) 
            case (neumann)
#if(_DIM_==1)            
               if (this%m%nf1min==0) then
                     this%uf(0) = f * this%uf(0)
               end if      
#elif(_DIM_==2)          
               if (this%m%nf1min==0) then
                     this%uf(0,:) = f * this%uf(0,:)
               end if      
               if (this%m%nf2min==0) then
                     this%uf(:,0) = f * this%uf(:,0)
               end if      
#elif(_DIM_==3)
               if (this%m%nf1min==0) then
                     this%uf(0,:,:) = f * this%uf(0,:,:)
               end if      
               if (this%m%nf2min==0) then
                     this%uf(:,0,:) = f * this%uf(:,0,:)
               end if      
               if (this%m%nf3min==0) then
                     this%uf(:,:,0) = f * this%uf(:,:,0)
               end if     
#endif
            end select

            select case(this%m%boundary_conditions) 
            case (periodic)
#ifdef _MPI_
#ifdef _REAL_
                call fftw_mpi_execute_dft_c2r(this%plan_backward, this%ucp, this%up)
#else
                call fftw_mpi_execute_dft(this%plan_backward, this%up, this%up)
#endif
#else
#ifdef _REAL_
                call fftw_execute_dft_c2r(this%plan_backward, this%ucp, this%up)
#else
                call fftw_execute_dft(this%plan_backward, this%up, this%up)
#endif
#endif

#ifndef _OPENMP
#if(_DIM_==1)
                this%u = this%u * (1.0_prec / this%m%g%nx)
#elif(_DIM_==2)
                this%u = this%u * (1.0_prec / (this%m%g%nx * this%m%g%ny))
#elif(_DIM_==3)
                this%u = this%u * (1.0_prec / (this%m%g%nx * this%m%g%ny * this%m%g%nz))
#endif

#else
!$OMP PARALLEL DO PRIVATE(j, u) 
                do j=1,n_threads
#if(_DIM_==1)
                   u => this%u(lbound(this%u,1)+this%m%g%jj(j-1):lbound(this%u,1)+this%m%g%jj(j)-1)
                   u = u * (1.0_prec / this%m%g%nx)
#elif(_DIM_==2)
                   u => this%u(:,lbound(this%u,2)+this%m%g%jj(j-1):lbound(this%u,2)+this%m%g%jj(j)-1)
                   u = u * (1.0_prec / (this%m%g%nx * this%m%g%ny))
#elif(_DIM_==3)
                   u => this%u(:,:,lbound(this%u,3)+this%m%g%jj(j-1):lbound(this%u,3)+this%m%g%jj(j)-1)
                   u = u * (1.0_prec / (this%m%g%nx * this%m%g%ny * this%m%g%nz))
#endif
                end do
!$OMP END PARALLEL DO 
#endif
            case (dirichlet, neumann)
#ifdef _MPI_
#ifdef _REAL_
                call fftw_mpi_execute_r2r(this%plan_backward, this%up, this%up)
#else
                call fftw_mpi_execute_r2r(this%plan_backward, this%urp, this%urp)
#endif
#else
#ifdef _REAL_
                call fftw_execute_r2r(this%plan_backward, this%up, this%up)
#else
                call fftw_execute_r2r(this%plan_backward, this%urp, this%urp)
#endif
#endif

#ifndef _OPENMP
#if(_DIM_==1)
                this%u = this%u * (1.0_prec / (2.0_prec*this%m%g%nx))
#elif(_DIM_==2)
                this%u = this%u * (1.0_prec / (4.0_prec * this%m%g%nx * this%m%g%ny))
#elif(_DIM_==3)
                this%u = this%u * (1.0_prec / (8.0_prec * this%m%g%nx * this%m%g%ny * this%m%g%nz))
#endif
#else
!$OMP PARALLEL DO PRIVATE(j, u) 
                do j=1,n_threads
#if(_DIM_==1)
                    u => this%u(lbound(this%u,1)+this%m%g%jj(j-1):lbound(this%u,1)+this%m%g%jj(j)-1)
                    u = u * (1.0_prec / (2.0_prec*this%m%g%nx))
#elif(_DIM_==2)
                    u => this%u(:,lbound(this%u,2)+this%m%g%jj(j-1):lbound(this%u,2)+this%m%g%jj(j)-1)
                    u = u * (1.0_prec / (4.0_prec * this%m%g%nx * this%m%g%ny))
#elif(_DIM_==3)
                    u => this%u(:,:,lbound(this%u,3)+this%m%g%jj(j-1):lbound(this%u,3)+this%m%g%jj(j)-1)
                    u = u * (1.0_prec / (8.0_prec * this%m%g%nx * this%m%g%ny * this%m%g%nz))
#endif
                end do
!$OMP END PARALLEL DO 
#endif
            end select

            this%is_real_space = .true.
        end if     
    end subroutine to_real_space



    subroutine propagate_A(this, dt)
        implicit none
        class(S(wf_fourier)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt

#ifdef _OPENMP
#if(_DIM_==1)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf(:)
        complex(kind=prec), pointer :: ufc(:)
#elif(_DIM_==2)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf(:,:)
        complex(kind=prec), pointer :: ufc(:,:)
#elif(_DIM_==3)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf(:,:,:)
        complex(kind=prec), pointer :: ufc(:,:,:)
#endif
        real(kind=prec), pointer :: ev(:)
        complex(kind=prec), pointer :: evc(:)
        integer :: j
#endif        

        call this%to_frequency_space

#if(_DIM_==1)
#ifdef _REAL_        
        select case(this%m%boundary_conditions) 
        case (periodic)
#ifndef _OPENMP
            this%ufc = exp((dt*this%coefficient)*this%m%eigenvalues1) * this%ufc
#else
!$OMP PARALLEL DO PRIVATE(j, ufc, ev) 
            do j=1,n_threads
                ufc => this%ufc(lbound(this%ufc,1)+this%m%jf(j-1):lbound(this%ufc,1)+this%m%jf(j)-1)
                ev => this%m%eigenvalues1(lbound(this%m%eigenvalues1,1)+this%m%jf(j-1):&
                                           lbound(this%m%eigenvalues1,1)+this%m%jf(j)-1)
                ufc = exp((dt*this%coefficient)*ev) * ufc
            end do
!$OMP END PARALLEL DO
#endif

        case (dirichlet, neumann)
#endif
#ifndef _OPENMP
            this%uf = exp((dt*this%coefficient)*this%m%eigenvalues1) * this%uf
#else 
!$OMP PARALLEL DO PRIVATE(j, uf, ev) 
            do j=1,n_threads
                uf => this%uf(lbound(this%uf,1)+this%m%jf(j-1):lbound(this%uf,1)+this%m%jf(j)-1)
                ev => this%m%eigenvalues1(lbound(this%m%eigenvalues1,1)+this%m%jf(j-1):&
                                           lbound(this%m%eigenvalues1,1)+this%m%jf(j)-1)
                uf = exp((dt*this%coefficient)*ev) * uf
            end do
!$OMP END PARALLEL DO 
#endif
#ifdef _REAL_        
        end select
#endif
#elif(_DIM_==2)
#ifdef _REAL_        
        select case(this%m%boundary_conditions) 
        case (periodic)
#ifndef _OPENMP
        this%ufc = spread(exp((dt*this%coefficient)*this%m%eigenvalues1),2, this%m%nf2max-this%m%nf2min+1) * this%ufc
        this%ufc = spread(exp((dt*this%coefficient)*this%m%eigenvalues2),1, this%m%nf1max-this%m%nf1min+1) * this%ufc
#else
!$OMP PARALLEL DO PRIVATE(j, ufc, ev) 
            do j=1,n_threads
                ufc => this%ufc(:,lbound(this%ufc,2)+this%m%jf(j-1):lbound(this%ufc,2)+this%m%jf(j)-1)
                ev => this%m%eigenvalues2(lbound(this%m%eigenvalues2,1)+this%m%jf(j-1):&
                                           lbound(this%m%eigenvalues2,1)+this%m%jf(j)-1)
                ufc = spread(exp((dt*this%coefficient)*this%m%eigenvalues1),2, this%m%jf(j)-this%m%jf(j-1)) * ufc
                ufc = spread(exp((dt*this%coefficient)*ev),1, this%m%nf1max-this%m%nf1min+1) * ufc
            end do
!$OMP END PARALLEL DO 
#endif
        case (dirichlet, neumann)
#endif

#ifndef _OPENMP
        this%uf = spread(exp((dt*this%coefficient)*this%m%eigenvalues1),2, this%m%nf2max-this%m%nf2min+1) * this%uf
        this%uf = spread(exp((dt*this%coefficient)*this%m%eigenvalues2),1, this%m%nf1max-this%m%nf1min+1) * this%uf
#else 
!$OMP PARALLEL DO PRIVATE(j, uf, ev) 
            do j=1,n_threads
                uf => this%uf(:,lbound(this%uf,2)+this%m%jf(j-1):lbound(this%uf,2)+this%m%jf(j)-1)
                ev => this%m%eigenvalues2(lbound(this%m%eigenvalues2,1)+this%m%jf(j-1):&
                                           lbound(this%m%eigenvalues2,1)+this%m%jf(j)-1)
                uf = spread(exp((dt*this%coefficient)*this%m%eigenvalues1),2, this%m%jf(j)-this%m%jf(j-1)) * uf
                uf = spread(exp((dt*this%coefficient)*ev),1, this%m%nf1max-this%m%nf1min+1) * uf
            end do
!$OMP END PARALLEL DO 
#endif
#ifdef _REAL_        
        end select
#endif
#elif(_DIM_==3)
#ifdef _REAL_        
        select case(this%m%boundary_conditions) 
        case (periodic)
#ifndef _OPENMP
        this%ufc = spread(spread(exp((dt*this%coefficient)*this%m%eigenvalues1), 2, this%m%nf2max-this%m%nf2min+1), &
                                3, this%m%nf3max-this%m%nf3min+1) * this%ufc
        this%ufc = spread(spread(exp((dt*this%coefficient)*this%m%eigenvalues2), 1, this%m%nf1max-this%m%nf1min+1), &
                                3, this%m%nf3max-this%m%nf3min+1) * this%ufc
        this%ufc = spread(spread(exp((dt*this%coefficient)*this%m%eigenvalues3), 1, this%m%nf1max-this%m%nf1min+1), &
                                2, this%m%nf2max-this%m%nf2min+1) * this%ufc
#else
!$OMP PARALLEL DO PRIVATE(j, ufc, ev) 
            do j=1,n_threads
                ufc => this%ufc(:,:,lbound(this%ufc,3)+this%m%jf(j-1):lbound(this%ufc,3)+this%m%jf(j)-1)
                ev => this%m%eigenvalues3(lbound(this%m%eigenvalues3,1)+this%m%jf(j-1):&
                                           lbound(this%m%eigenvalues3,1)+this%m%jf(j)-1)
                ufc = spread(spread(exp((dt*this%coefficient)*this%m%eigenvalues1), &
                                2, this%m%nf2max-this%m%nf2min+1), &
                                3,  this%m%jf(j)-this%m%jf(j-1)) * ufc
                ufc = spread(spread(exp((dt*this%coefficient)*this%m%eigenvalues2), &
                                1, this%m%nf1max-this%m%nf1min+1), &
                                3,  this%m%jf(j)-this%m%jf(j-1)) * ufc
                ufc = spread(spread(exp((dt*this%coefficient)*ev), &
                                1, this%m%nf1max-this%m%nf1min+1), &
                                2, this%m%nf2max-this%m%nf2min+1) * ufc
            end do
!$OMP END PARALLEL DO 
#endif
        case (dirichlet, neumann)
#endif
#ifndef _OPENMP
        this%uf = spread(spread(exp((dt*this%coefficient)*this%m%eigenvalues1), 2, this%m%nf2max-this%m%nf2min+1), &
                                3, this%m%nf3max-this%m%nf3min+1) * this%uf
        this%uf = spread(spread(exp((dt*this%coefficient)*this%m%eigenvalues2), 1, this%m%nf1max-this%m%nf1min+1), &
                                3, this%m%nf3max-this%m%nf3min+1) * this%uf
        this%uf = spread(spread(exp((dt*this%coefficient)*this%m%eigenvalues3), 1, this%m%nf1max-this%m%nf1min+1), &
                                2, this%m%nf2max-this%m%nf2min+1) * this%uf
#else
!$OMP PARALLEL DO PRIVATE(j, uf, ev) 
            do j=1,n_threads
                uf => this%uf(:,:,lbound(this%uf,3)+this%m%jf(j-1):lbound(this%uf,3)+this%m%jf(j)-1)
                ev => this%m%eigenvalues3(lbound(this%m%eigenvalues3,1)+this%m%jf(j-1):&
                                           lbound(this%m%eigenvalues3,1)+this%m%jf(j)-1)
                uf = spread(spread(exp((dt*this%coefficient)*this%m%eigenvalues1), &
                                2, this%m%nf2max-this%m%nf2min+1), &
                                3,  this%m%jf(j)-this%m%jf(j-1)) * uf
                uf = spread(spread(exp((dt*this%coefficient)*this%m%eigenvalues2), &
                                1, this%m%nf1max-this%m%nf1min+1), &
                                3,  this%m%jf(j)-this%m%jf(j-1)) * uf
                uf = spread(spread(exp((dt*this%coefficient)*ev), &
                                1, this%m%nf1max-this%m%nf1min+1), &
                                2, this%m%nf2max-this%m%nf2min+1) * uf
            end do
!$OMP END PARALLEL DO 
#endif
#ifdef _REAL_        
        end select
#endif
#endif
    end subroutine propagate_A



    subroutine add_apply_A(this, wf, coefficient)
        class(S(wf_fourier)), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: wf 
        _COMPLEX_OR_REAL_(kind=prec), intent(in), optional :: coefficient
        _COMPLEX_OR_REAL_(kind=prec) :: C 

#ifdef _OPENMP
#if(_DIM_==1)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf1(:)
        complex(kind=prec), pointer :: ufc1(:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf2(:)
        complex(kind=prec), pointer :: ufc2(:)
#elif(_DIM_==2)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf1(:,:)
        complex(kind=prec), pointer :: ufc1(:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf2(:,:)
        complex(kind=prec), pointer :: ufc2(:,:)
#elif(_DIM_==3)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf1(:,:,:)
        complex(kind=prec), pointer :: ufc1(:,:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf2(:,:,:)
        complex(kind=prec), pointer :: ufc2(:,:,:)
#endif
        real(kind=prec), pointer :: ev(:)
        complex(kind=prec), pointer :: evc(:)
        integer :: j
#endif        
        

        select type (wf)
        class is (S(wf_fourier))
        
        if (.not.associated(wf%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    

        C = this%coefficient
        if (present(coefficient)) then
            C = C*coefficient
        end if    

        call this%to_frequency_space
        call wf%to_frequency_space

#if(_DIM_==1)
#ifdef _REAL_        
        select case(this%m%boundary_conditions) 
        case (periodic)
#ifndef _OPENMP            
            wf%ufc = wf%ufc + C*this%m%eigenvalues1*this%ufc
#else
!$OMP PARALLEL DO PRIVATE(j, ufc1, ufc2, ev) 
            do j=1,n_threads
                ufc1 => wf%ufc(lbound(wf%ufc,1)+this%m%jf(j-1):lbound(wf%ufc,1)+this%m%jf(j)-1)
                ufc2 => this%ufc(lbound(this%ufc,1)+this%m%jf(j-1):lbound(this%ufc,1)+this%m%jf(j)-1)
                ev => this%m%eigenvalues1(lbound(this%m%eigenvalues1,1)+this%m%jf(j-1):&
                                           lbound(this%m%eigenvalues1,1)+this%m%jf(j)-1)
                ufc1 = ufc1 + C*ev*ufc2
            end do
!$OMP END PARALLEL DO
#endif
        case (dirichlet, neumann)
#endif
#ifndef _OPENMP            
            wf%uf = wf%uf + C*this%m%eigenvalues1*this%uf
#else
!$OMP PARALLEL DO PRIVATE(j, uf1, uf2, ev) 
            do j=1,n_threads
                uf1 => wf%uf(lbound(wf%uf,1)+this%m%jf(j-1):lbound(wf%uf,1)+this%m%jf(j)-1)
                uf2 => this%uf(lbound(this%uf,1)+this%m%jf(j-1):lbound(this%uf,1)+this%m%jf(j)-1)
                ev => this%m%eigenvalues1(lbound(this%m%eigenvalues1,1)+this%m%jf(j-1):&
                                           lbound(this%m%eigenvalues1,1)+this%m%jf(j)-1)
                uf1 = uf1 + C*ev*uf2
            end do
!$OMP END PARALLEL DO 
#endif
#ifdef _REAL_        
        end select
#endif
#elif(_DIM_==2)
#ifdef _REAL_        
        select case(this%m%boundary_conditions) 
        case (periodic)
#ifndef _OPENMP            
        wf%ufc = wf%ufc + C*spread(this%m%eigenvalues1, 2, this%m%nf2max-this%m%nf2min+1) * this%ufc
        wf%ufc = wf%ufc + C*spread(this%m%eigenvalues2, 1, this%m%nf1max-this%m%nf1min+1) * this%ufc
#else
!$OMP PARALLEL DO PRIVATE(j, ufc1, ufc2, ev) 
            do j=1,n_threads
                ufc1 => wf%ufc(:,lbound(wf%ufc,2)+this%m%jf(j-1):lbound(wf%ufc,2)+this%m%jf(j)-1)
                ufc2 => this%ufc(:,lbound(this%ufc,2)+this%m%jf(j-1):lbound(this%ufc,2)+this%m%jf(j)-1)
                ev => this%m%eigenvalues2(lbound(this%m%eigenvalues2,1)+this%m%jf(j-1):&
                                           lbound(this%m%eigenvalues2,1)+this%m%jf(j)-1)
                ufc1 = ufc1 + C*spread(this%m%eigenvalues1, 2,  this%m%jf(j)-this%m%jf(j-1)) * ufc2
                ufc1 = ufc1 + C*spread(ev, 1, this%m%nf1max-this%m%nf1min+1) * ufc2
            end do
!$OMP END PARALLEL DO 
#endif
        case (dirichlet, neumann)
#endif
#ifndef _OPENMP            
        wf%uf = wf%uf + C*spread(this%m%eigenvalues1,2, this%m%nf2max-this%m%nf2min+1) * this%uf
        wf%uf = wf%uf + C*spread(this%m%eigenvalues2,1, this%m%nf1max-this%m%nf1min+1) * this%uf
#else
!$OMP PARALLEL DO PRIVATE(j, uf1, uf2, ev) 
            do j=1,n_threads
                uf1 => wf%uf(:,lbound(wf%uf,2)+this%m%jf(j-1):lbound(wf%uf,2)+this%m%jf(j)-1)
                uf2 => this%uf(:,lbound(this%uf,2)+this%m%jf(j-1):lbound(this%uf,2)+this%m%jf(j)-1)
                ev => this%m%eigenvalues2(lbound(this%m%eigenvalues2,1)+this%m%jf(j-1):&
                                           lbound(this%m%eigenvalues2,1)+this%m%jf(j)-1)
                uf1 = uf1 + C*spread(this%m%eigenvalues1, 2,  this%m%jf(j)-this%m%jf(j-1)) * uf2
                uf1 = uf1 + C*spread(ev, 1, this%m%nf1max-this%m%nf1min+1) * uf2
            end do
!$OMP END PARALLEL DO 
#endif
#ifdef _REAL_        
        end select
#endif
#elif(_DIM_==3)
#ifdef _REAL_        
        select case(this%m%boundary_conditions) 
        case (periodic)
#ifndef _OPENMP            
        wf%ufc = wf%ufc + C*spread(spread(this%m%eigenvalues1, 2, this%m%nf2max-this%m%nf2min+1), &
                                3, this%m%nf3max-this%m%nf3min+1) * this%ufc
        wf%ufc = wf%ufc + C*spread(spread(this%m%eigenvalues2, 1, this%m%nf1max-this%m%nf1min+1), &
                                3, this%m%nf3max-this%m%nf3min+1) * this%ufc
        wf%ufc = wf%ufc + C*spread(spread(this%m%eigenvalues3, 1, this%m%nf1max-this%m%nf1min+1), &
                                2, this%m%nf2max-this%m%nf2min+1) * this%ufc
#else
!$OMP PARALLEL DO PRIVATE(j, ufc1, ufc2, ev) 
            do j=1,n_threads
                ufc1 => wf%ufc(:,:,lbound(wf%ufc,3)+this%m%jf(j-1):lbound(wf%ufc,3)+this%m%jf(j)-1)
                ufc2 => this%ufc(:,:,lbound(this%ufc,3)+this%m%jf(j-1):lbound(this%ufc,3)+this%m%jf(j)-1)
                ev => this%m%eigenvalues3(lbound(this%m%eigenvalues3,1)+this%m%jf(j-1):&
                                           lbound(this%m%eigenvalues3,1)+this%m%jf(j)-1)
                ufc1 = ufc1 + C*spread(spread(this%m%eigenvalues1, &
                                2, this%m%nf2max-this%m%nf2min+1), &
                                3, this%m%jf(j)-this%m%jf(j-1)) * ufc2
                ufc1 = ufc1 + C*spread(spread(this%m%eigenvalues2, &
                                1, this%m%nf1max-this%m%nf1min+1), &
                                3, this%m%jf(j)-this%m%jf(j-1)) * ufc2
                ufc1 = ufc1 + C*spread(spread(ev, &
                                1, this%m%nf1max-this%m%nf1min+1), &
                                2, this%m%nf2max-this%m%nf2min+1) * ufc2
            end do
!$OMP END PARALLEL DO 
#endif
        case (dirichlet, neumann)
#endif
#ifndef _OPENMP  
        wf%uf = wf%uf + C*spread(spread(this%m%eigenvalues1, 2, this%m%nf2max-this%m%nf2min+1), &
                                3, this%m%nf3max-this%m%nf3min+1) * this%uf
        wf%uf = wf%uf + C*spread(spread(this%m%eigenvalues2, 1, this%m%nf1max-this%m%nf1min+1), &
                                3, this%m%nf3max-this%m%nf3min+1) * this%uf
        wf%uf = wf%uf + C*spread(spread(this%m%eigenvalues3, 1, this%m%nf1max-this%m%nf1min+1), &
                                2, this%m%nf2max-this%m%nf2min+1) * this%uf
#else
!$OMP PARALLEL DO PRIVATE(j, uf1, uf2, ev) 
            do j=1,n_threads
                uf1 => wf%uf(:,:,lbound(wf%uf,3)+this%m%jf(j-1):lbound(wf%uf,3)+this%m%jf(j)-1)
                uf2 => this%uf(:,:,lbound(this%uf,3)+this%m%jf(j-1):lbound(this%uf,3)+this%m%jf(j)-1)
                ev => this%m%eigenvalues3(lbound(this%m%eigenvalues3,1)+this%m%jf(j-1):&
                                           lbound(this%m%eigenvalues3,1)+this%m%jf(j)-1)
                uf1 = uf1 + C*spread(spread(this%m%eigenvalues1, &
                                2, this%m%nf2max-this%m%nf2min+1), &
                                3, this%m%jf(j)-this%m%jf(j-1)) * uf2
                uf1 = uf1 + C*spread(spread(this%m%eigenvalues2, &
                                1, this%m%nf1max-this%m%nf1min+1), &
                                3, this%m%jf(j)-this%m%jf(j-1)) * uf2
                uf1 = uf1 + C*spread(spread(ev, &
                                1, this%m%nf1max-this%m%nf1min+1), &
                                2, this%m%nf2max-this%m%nf2min+1) * uf2

            end do
!$OMP END PARALLEL DO 
#endif
#ifdef _REAL_        
        end select
#endif
#endif
        end select
    end subroutine add_apply_A

   
    subroutine save(this, filename)
#ifdef _NO_HDF5_
        class(S(wf_fourier)), intent(inout) :: this
        character(len=*), intent(in) :: filename
        print *, "W: save not implemented"
#else
#ifdef _QUADPRECISION_
        use tssmq_hdf5
#else
        use tssm_hdf5
#endif        
        class(S(wf_fourier)), intent(inout) :: this
        character(len=*), intent(in) :: filename
        
        call this%to_real_space

#ifdef _REAL_        
        call hdf5_save_real_gridfun(this%m%g, this%u, filename, trim(this%m%dset_name))
#else        
        call hdf5_save_complex_gridfun(this%m%g, this%u, filename, &
                     trim(this%m%dset_name_real), trim(this%m%dset_name_imag))
#endif        
        call hdf5_write_grid_attributes(this%m%g, filename)
#endif        
    end subroutine save



    subroutine load(this, filename)
#ifdef _NO_HDF5_
        class(S(wf_fourier)), intent(inout) :: this
        character(len=*), intent(in) :: filename
        print *, "W: load not implemented"
#else
#ifdef _QUADPRECISION_
        use tssmq_hdf5
#else
        use tssm_hdf5
#endif        
        class(S(wf_fourier)), intent(inout) :: this
        character(len=*), intent(in) :: filename
#ifdef _QUADPRECISION_
        real(kind=prec), parameter :: eps = epsilon(1.0_8)
#else        
        real(kind=prec), parameter :: eps = epsilon(1.0_prec)
#endif
        integer :: offset(3)
#if(_DIM_==1)
        type(grid_equidistant_1D) :: g 
 
#elif(_DIM_==2)
        type(grid_equidistant_2D) :: g 
#elif(_DIM_==3)
        type(grid_equidistant_3D) :: g 
#endif 
#ifdef _REAL_
        call hdf5_read_grid_attributes(g, filename, trim(this%m%dset_name))
#else        
        call hdf5_read_grid_attributes(g, filename, trim(this%m%dset_name_real))
#endif        
#if(_DIM_==1)
        if (.not.(abs(this%m%g%dx-g%dx)<eps)) then
            stop "E: incompatible grids (steps sizes do not match)"
        end if
        offset(1) = nint((g%xmin-this%m%g%xmin)/this%m%g%dx)
        this%m%g%xmin=g%xmin-offset(1)*this%m%g%dx
        this%m%g%xmax=this%m%g%xmin+this%m%g%nx*this%m%g%dx
#elif(_DIM_==2)
        if (.not.(abs(this%m%g%dx-g%dx)<eps.and.abs(this%m%g%dy-g%dy)<eps)) then
            stop "E: incompatible grids (steps sizes do not match)"
        end if
        offset(1) = nint((g%xmin-this%m%g%xmin)/this%m%g%dx)
        this%m%g%xmin=g%xmin-offset(1)*this%m%g%dx
        this%m%g%xmax=this%m%g%xmin+this%m%g%nx*this%m%g%dx
        offset(2) = nint((g%ymin-this%m%g%ymin)/this%m%g%dy)
        this%m%g%ymin=g%ymin-offset(2)*this%m%g%dy
        this%m%g%ymax=this%m%g%ymin+this%m%g%ny*this%m%g%dy
#elif(_DIM_==3)
        if (.not.(abs(this%m%g%dx-g%dx)<eps.and.abs(this%m%g%dy-g%dy)<eps.and.abs(this%m%g%dz-g%dz)<eps)) then
            stop "E: incompatible grids (steps sizes do not match)"
        end if
        offset(1) = nint((g%xmin-this%m%g%xmin)/this%m%g%dx)
        this%m%g%xmin=g%xmin-offset(1)*this%m%g%dx
        this%m%g%xmax=this%m%g%xmin+this%m%g%nx*this%m%g%dx
        offset(2) = nint((g%ymin-this%m%g%ymin)/this%m%g%dy)
        this%m%g%ymin=g%ymin-offset(2)*this%m%g%dy
        this%m%g%ymax=this%m%g%ymin+this%m%g%ny*this%m%g%dy
        offset(3) = nint((g%zmin-this%m%g%zmin)/this%m%g%dz)
        this%m%g%zmin=g%zmin-offset(3)*this%m%g%dz
        this%m%g%zmax=this%m%g%zmin+this%m%g%nz*this%m%g%dz
#endif
        this%u = 0.0_prec
#ifdef _REAL_        
        call hdf5_load_real_gridfun(this%m%g, this%u, filename, trim(this%m%dset_name), offset=offset)
#else      
        call hdf5_load_complex_gridfun(this%m%g, this%u, filename, & 
                     trim(this%m%dset_name_real), trim(this%m%dset_name_imag), offset=offset)
#endif        
        this%is_real_space = .true.
#endif
    end subroutine load


    subroutine set(this, f)
        class(S(wf_fourier)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), external :: f
#ifdef _REAL_        
        call this%m%g%set_real_gridfun(this%u, f)
#else
        call this%m%g%set_complex_gridfun(this%u, f)
#endif        
        this%is_real_space = .true.
    end subroutine set

    subroutine set_t(this, f, t)
        class(S(wf_fourier)), intent(inout) :: this
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
        class(S(wf_fourier)), intent(inout) :: this
        real(kind=prec), external :: f
        call this%m%g%rset_complex_gridfun(this%u, f)
        this%is_real_space = .true.
    end subroutine rset

    subroutine rset_t(this, f, t)
        class(S(wf_fourier)), intent(inout) :: this
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t
        call this%m%g%rset_t_complex_gridfun(this%u, f, t)
        this%is_real_space = .true.
    end subroutine rset_t

#endif        

    function norm(this) result(n)
        class(S(wf_fourier)), intent(inout) :: this
        real(kind=prec) :: n
        
        !!! TODO handle norm in frequency space without transforming
        call this%to_real_space 
#ifdef _REAL_        
        n = this%m%g%norm_real_gridfun(this%u)
#else
        n = this%m%g%norm_complex_gridfun(this%u)
#endif        
    end function norm


    subroutine normalize(this, norm)
        class(S(wf_fourier)), intent(inout) :: this
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



    subroutine scale(this, factor)
        class(S(wf_fourier)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: factor
#ifdef _OPENMP
#if(_DIM_==1)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:)
#elif(_DIM_==2)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:)
#elif(_DIM_==3)            
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
#if(_DIM_==1)
              u => this%u(lbound(this%u,1)+this%m%g%jj(j-1):lbound(this%u,1)+this%m%g%jj(j)-1)
#elif(_DIM_==2)
              u => this%u(:,lbound(this%u,2)+this%m%g%jj(j-1):lbound(this%u,2)+this%m%g%jj(j)-1)
#elif(_DIM_==3)
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


    subroutine axpy(this, other, factor)
         class(S(wf_fourier)), intent(inout) :: this
         class(_WAVE_FUNCTION_), intent(inout) :: other
         _COMPLEX_OR_REAL_(kind=prec), intent(in) :: factor
#ifdef _OPENMP
#if(_DIM_==1)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u1(:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u2(:)
#elif(_DIM_==2)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u1(:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u2(:,:)
#elif(_DIM_==3)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u1(:,:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u2(:,:,:)
#endif
        integer :: j
#endif        
        select type (other)
        class is (S(wf_fourier))
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
#if(_DIM_==1)
              u1 => this%u(lbound(this%u,1)+this%m%g%jj(j-1):lbound(this%u,1)+this%m%g%jj(j)-1)
              u2 => other%u(lbound(other%u,1)+other%m%g%jj(j-1):lbound(other%u,1)+other%m%g%jj(j)-1)
#elif(_DIM_==2)
              u1 => this%u(:,lbound(this%u,2)+this%m%g%jj(j-1):lbound(this%u,2)+this%m%g%jj(j)-1)
              u2 => other%u(:,lbound(other%u,2)+other%m%g%jj(j-1):lbound(other%u,2)+other%m%g%jj(j)-1)
#elif(_DIM_==3)
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



    subroutine copy(this, source)
        class(S(wf_fourier)), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: source 
#ifdef _OPENMP
#if(_DIM_==1)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u1(:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u2(:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf1(:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf2(:)
        complex(kind=prec), pointer :: ufc1(:)
        complex(kind=prec), pointer :: ufc2(:)
#elif(_DIM_==2)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u1(:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u2(:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf1(:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf2(:,:)
        complex(kind=prec), pointer :: ufc1(:,:)
        complex(kind=prec), pointer :: ufc2(:,:)
#elif(_DIM_==3)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u1(:,:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u2(:,:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf1(:,:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf2(:,:,:)
        complex(kind=prec), pointer :: ufc1(:,:,:)
        complex(kind=prec), pointer :: ufc2(:,:,:)
#endif
        integer :: j
#endif        
        select type (source)
        class is (S(wf_fourier))
        if (.not.associated(source%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    
       
        this%is_real_space = source%is_real_space

        if (this%is_real_space) then
#ifndef _OPENMP
            this%u = source%u
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
             do j=1,n_threads
#if(_DIM_==1)
              u1 => this%u(lbound(this%u,1)+this%m%g%jj(j-1):lbound(this%u,1)+this%m%g%jj(j)-1)
              u2 => source%u(lbound(source%u,1)+this%m%g%jj(j-1):lbound(source%u,1)+this%m%g%jj(j)-1)
#elif(_DIM_==2)
              u1 => this%u(:,lbound(this%u,2)+this%m%g%jj(j-1):lbound(this%u,2)+this%m%g%jj(j)-1)
              u2 => source%u(:,lbound(source%u,2)+this%m%g%jj(j-1):lbound(source%u,2)+this%m%g%jj(j)-1)
#elif(_DIM_==3)
              u1 => this%u(:,:,lbound(this%u,3)+this%m%g%jj(j-1):lbound(this%u,3)+this%m%g%jj(j)-1)
              u2 => source%u(:,:,lbound(source%u,3)+source%m%g%jj(j-1):lbound(source%u,3)+this%m%g%jj(j)-1)
#endif
              u1 = u2
             end do
!$OMP END PARALLEL DO 
#endif
        else    
#ifdef _REAL_        
            select case(this%m%boundary_conditions) 
            case (periodic)
#ifndef _OPENMP
                this%ufc = source%ufc
#else
!$OMP PARALLEL DO PRIVATE(j, ufc1, ufc2) 
             do j=1,n_threads
#if(_DIM_==1)
              ufc1 => this%ufc(lbound(this%ufc,1)+this%m%jf(j-1):lbound(this%ufc,1)+this%m%jf(j)-1)
              ufc2 => source%ufc(lbound(source%ufc,1)+this%m%jf(j-1):lbound(source%ufc,1)+this%m%jf(j)-1)
#elif(_DIM_==2)
              ufc1 => this%ufc(:,lbound(this%ufc,2)+this%m%jf(j-1):lbound(this%ufc,2)+this%m%jf(j)-1)
              ufc2 => source%ufc(:,lbound(source%ufc,2)+this%m%jf(j-1):lbound(source%ufc,2)+this%m%jf(j)-1)
#elif(_DIM_==3)
              ufc1 => this%ufc(:,:,lbound(this%ufc,3)+this%m%jf(j-1):lbound(this%ufc,3)+this%m%jf(j)-1)
              ufc2 => source%ufc(:,:,lbound(source%ufc,3)+source%m%jf(j-1):lbound(source%ufc,3)+this%m%jf(j)-1)
#endif
              ufc1 = ufc2
             end do
!$OMP END PARALLEL DO 
#endif
            case (dirichlet, neumann)
#endif
!TODO OMP
#ifndef _OPENMP
                this%uf = source%uf 
#else
!$OMP PARALLEL DO PRIVATE(j, uf1, uf2) 
             do j=1,n_threads
#if(_DIM_==1)
              uf1 => this%uf(lbound(this%uf,1)+this%m%jf(j-1):lbound(this%uf,1)+this%m%jf(j)-1)
              uf2 => source%uf(lbound(source%uf,1)+this%m%jf(j-1):lbound(source%uf,1)+this%m%jf(j)-1)
#elif(_DIM_==2)
              uf1 => this%uf(:,lbound(this%uf,2)+this%m%jf(j-1):lbound(this%uf,2)+this%m%jf(j)-1)
              uf2 => source%uf(:,lbound(source%uf,2)+this%m%jf(j-1):lbound(source%uf,2)+this%m%jf(j)-1)
#elif(_DIM_==3)
              uf1 => this%uf(:,:,lbound(this%uf,3)+this%m%jf(j-1):lbound(this%uf,3)+this%m%jf(j)-1)
              uf2 => source%uf(:,:,lbound(source%uf,3)+source%m%jf(j-1):lbound(source%uf,3)+this%m%jf(j)-1)
#endif
              uf1 = uf2
             end do
!$OMP END PARALLEL DO 
#endif
#ifdef _REAL_        
            end select
#endif
        end if

        class default
           stop "E: wave functions not belonging to the same method"
        end select
    end subroutine copy


    function distance(this, wf) result(n)
        class(S(wf_fourier)), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: wf 
        real(kind=prec) :: n
#ifdef _OPENMP
#if(_DIM_==1)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf1(:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf2(:)
#elif(_DIM_==2)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf1(:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf2(:,:)
#elif(_DIM_==3)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf1(:,:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf2(:,:,:)
#endif
        integer :: j
#endif        
        select type (wf)
        class is (S(wf_fourier))
        if (.not.associated(wf%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    
        if (associated(this%u,wf%u)) then
             n = 0.0_prec
             return
        end if

        !!! TODO handle norm in frequency space without transforming
        call this%to_frequency_space 
        call wf%to_frequency_space 
#ifndef _OPENMP        
        this%uf = this%uf - wf%uf
#else
!$OMP PARALLEL DO PRIVATE(j, uf1, uf2) 
        do j=1,n_threads
#if(_DIM_==1)
            uf1 => wf%uf(lbound(wf%uf,1)+this%m%jf(j-1):lbound(wf%uf,1)+this%m%jf(j)-1)
            uf2 => this%uf(lbound(this%uf,1)+this%m%jf(j-1):lbound(this%uf,1)+this%m%jf(j)-1)
#elif(_DIM_==2)
            uf1 => wf%uf(:,lbound(wf%uf,2)+this%m%jf(j-1):lbound(wf%uf,2)+this%m%jf(j)-1)
            uf2 => this%uf(:,lbound(this%uf,2)+this%m%jf(j-1):lbound(this%uf,2)+this%m%jf(j)-1)
#elif(_DIM_==3)
            uf1 => wf%uf(:,:,lbound(wf%uf,3)+this%m%jf(j-1):lbound(wf%uf,3)+this%m%jf(j)-1)
            uf2 => this%uf(:,:,lbound(this%uf,3)+this%m%jf(j-1):lbound(this%uf,3)+this%m%jf(j)-1)
#endif
            uf2 = uf2 - uf1
        end do
!$OMP END PARALLEL DO 
#endif
        n = this%norm_in_frequency_space()
#ifndef _OPENMP        
        this%uf = this%uf + wf%uf
#else
!$OMP PARALLEL DO PRIVATE(j, uf1, uf2) 
        do j=1,n_threads
#if(_DIM_==1)
            uf1 => wf%uf(lbound(wf%uf,1)+this%m%jf(j-1):lbound(wf%uf,1)+this%m%jf(j)-1)
            uf2 => this%uf(lbound(this%uf,1)+this%m%jf(j-1):lbound(this%uf,1)+this%m%jf(j)-1)
#elif(_DIM_==2)
            uf1 => wf%uf(:,lbound(wf%uf,2)+this%m%jf(j-1):lbound(wf%uf,2)+this%m%jf(j)-1)
            uf2 => this%uf(:,lbound(this%uf,2)+this%m%jf(j-1):lbound(this%uf,2)+this%m%jf(j)-1)
#elif(_DIM_==3)
            uf1 => wf%uf(:,:,lbound(wf%uf,3)+this%m%jf(j-1):lbound(wf%uf,3)+this%m%jf(j)-1)
            uf2 => this%uf(:,:,lbound(this%uf,3)+this%m%jf(j-1):lbound(this%uf,3)+this%m%jf(j)-1)
#endif
            uf2 = uf2 + uf1
        end do
!$OMP END PARALLEL DO 
#endif        
        class default
           stop "E: wave functions not belonging to the same method"
        end select
    end function distance



    function norm_in_frequency_space(this) result(N)
        class(S(wf_fourier)), intent(inout) :: this
        real(kind=prec) :: N

        real(kind=prec) :: h 
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: h1 
#endif
#ifdef _OPENMP
#if(_DIM_==1)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf(:)
        complex(kind=prec), pointer :: ufc(:)
#elif(_DIM_==2)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf(:,:)
        complex(kind=prec), pointer :: ufc(:,:)
#elif(_DIM_==3)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf(:,:,:)
        complex(kind=prec), pointer :: ufc(:,:,:)
#endif
        integer :: j
#endif        

        
        call this%to_frequency_space 
#ifdef _REAL_
        select case(this%m%boundary_conditions)
        case(periodic) !TODO correct for the factor 2
#ifndef _OPENMP
            h = 2.0_prec*sum(real(this%ufc,prec)**2 + aimag(this%ufc)**2) 
#else
            h = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, ufc) REDUCTION(+:h) 
            do j=1,n_threads
#if(_DIM_==1)
                ufc => this%ufc(lbound(this%ufc,1)+this%m%jf(j-1):lbound(this%ufc,1)+this%m%jf(j)-1)
#elif(_DIM_==2)
                ufc => this%ufc(:,lbound(this%ufc,2)+this%m%jf(j-1):lbound(this%ufc,2)+this%m%jf(j)-1)
#elif(_DIM_==3)
                ufc => this%ufc(:,:,lbound(this%ufc,3)+this%m%jf(j-1):lbound(this%ufc,3)+this%m%jf(j)-1)
#endif
                h = h + 2.0_prec*sum(real(ufc,prec)**2 + aimag(ufc)**2)
        end do
!$OMP END PARALLEL DO 
#endif
        case(dirichlet, neumann)
#ifndef _OPENMP
            h = sum( this%uf**2 )
#else
            h = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, uf) REDUCTION(+:h) 
            do j=1,n_threads
#if(_DIM_==1)
                uf => this%uf(lbound(this%uf,1)+this%m%jf(j-1):lbound(this%uf,1)+this%m%jf(j)-1)
#elif(_DIM_==2)
                uf => this%uf(:,lbound(this%uf,2)+this%m%jf(j-1):lbound(this%uf,2)+this%m%jf(j)-1)
#elif(_DIM_==3)
                uf => this%uf(:,:,lbound(this%uf,3)+this%m%jf(j-1):lbound(this%uf,3)+this%m%jf(j)-1)
#endif
                h = h + sum( uf**2 )
        end do
!$OMP END PARALLEL DO 
            
#endif
        end select
#else
#ifndef _OPENMP
        h = sum(real(this%uf,prec)**2 + aimag(this%uf)**2) 
#else
        h = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, uf) REDUCTION(+:h) 
        do j=1,n_threads
#if(_DIM_==1)
            uf => this%uf(lbound(this%uf,1)+this%m%jf(j-1):lbound(this%uf,1)+this%m%jf(j)-1)
#elif(_DIM_==2)
            uf => this%uf(:,lbound(this%uf,2)+this%m%jf(j-1):lbound(this%uf,2)+this%m%jf(j)-1)
#elif(_DIM_==3)
            uf => this%uf(:,:,lbound(this%uf,3)+this%m%jf(j-1):lbound(this%uf,3)+this%m%jf(j)-1)
#endif
            h = h + sum(real(uf,prec)**2 + aimag(uf)**2)
        end do
!$OMP END PARALLEL DO 
        
#endif
#endif
#ifdef _MPI_
        call MPI_Reduce(h, h1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            h = h1
#endif
#if(_DIM_==1)
            N =  h * this%m%g%dx/(this%m%g%nx)                                        
#elif(_DIM_==2)
            N =  h * this%m%g%dx*this%m%g%dy/(this%m%g%nx*this%m%g%ny)                                         
#elif(_DIM_==3)
            N =  h * this%m%g%dx*this%m%g%dy*this%m%g%dz/(this%m%g%nx*this%m%g%ny*this%m%g%nz) 
#endif
            select case(this%m%boundary_conditions)
            case(periodic)
            case(dirichlet, neumann)
                N = N/(2**_DIM_)
            end select
#ifdef _MPI_
        end if 
        call MPI_Bcast(N, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
!TODO check for sqrt !!!
        N = sqrt(N)
    end function norm_in_frequency_space

#ifdef _QUADPRECISION_
end module S(tssmq_fourier)
#else
end module S(tssm_fourier)
#endif


