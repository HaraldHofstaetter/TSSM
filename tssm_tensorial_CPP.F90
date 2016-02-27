!!! TODO: OPENMP- Version of tssm_tensorial_CPP.F90 (make "!$OMP PARALLEL WORKSHARE" explicit)
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
#define fftw_free fftwq_free
#endif

#ifdef _QUADPRECISION_
module S(tssmq_tensorial)
    use tssmq_base
    use tssmq_grid
    use tssmq_fourier_common
#else
module S(tssm_tensorial)
    use tssm_base
    use tssm_grid
    use tssm_fourier_common
#endif    
    implicit none

    private
    public ::  S(tensorial),  S(wf_tensorial)

    type, extends(spectral_method) :: S(tensorial)
#if(_DIM_==1)
        type(grid_tensorial_1D) :: g 
#elif(_DIM_==2)
        type(grid_tensorial_2D) :: g 
#elif(_DIM_==3)
        type(grid_tensorial_3D) :: g 
#endif 
        integer :: nf1
        integer :: nf1min !integer range for transformed data
        integer :: nf1max
#if(_DIM_>=2)
        integer :: nf2
        integer :: nf2min
        integer :: nf2max
#endif
#if(_DIM_>=3)
        integer :: nf3
        integer :: nf3min
        integer :: nf3max
#endif
        real(kind=prec), allocatable :: eigenvalues1(:)
        real(kind=prec), allocatable :: H1(:,:)
#if(_DIM_>=2)
        real(kind=prec), allocatable :: eigenvalues2(:)
        real(kind=prec), allocatable :: H2(:,:)
#endif 
#if(_DIM_>=3)
        real(kind=prec), allocatable :: eigenvalues3(:)
        real(kind=prec), allocatable :: H3(:,:)
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
        logical :: propagate_time_together_with_A = .true.
    contains    
        procedure :: finalize => finalize_method
        !final :: final_tensorial_1D
        !! Fortran 2003 feature final seems to be not properly implemented
        !! in the gcc/gfortran compiler :(
    end type S(tensorial)

    interface  S(tensorial) ! constructor
        module procedure new_method
    end interface S(tensorial)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    type, extends(_WAVE_FUNCTION_) :: S(wf_tensorial)
        class(S(tensorial)), pointer :: m 

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

        _COMPLEX_OR_REAL_(kind=prec) :: coefficient = 1.0_prec

    contains
        procedure :: to_real_space
        procedure :: to_frequency_space
        procedure :: propagate_A
        procedure :: propagate_A_derivative
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
    end type S(wf_tensorial)

    interface S(wf_tensorial) ! constructor
        module procedure new_wf
    end interface S(wf_tensorial)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains



#if(_DIM_==1)
    function new_method(n1 ) result(this)
#elif(_DIM_==2)
    function new_method(n1, n2) result(this)
#elif(_DIM_==3)
    function new_method(n1, n2, n3) result(this)
#endif
        type(S(tensorial)) :: this
        integer, intent(in) :: n1
#if(_DIM_>=2)
        integer, intent(in) :: n2
#endif
#if(_DIM_>=3)
        integer, intent(in) :: n3
#endif

        integer :: i

!adjust grid parameters, TODO MPI!!!


        this%g%nx = n1
        this%g%nn1min = 0
        this%g%nn1max = this%g%nx - 1

        this%g%n1min = this%g%nn1min
        this%g%n1max = this%g%nn1max
        this%nf1 = this%g%nx
        this%nf1min = this%g%nn1min
        this%nf1max = this%g%nn1max
        this%g%m1min = this%g%nn1min
        this%g%m1max = this%g%nn1max
#if(_DIM_>=2)
        this%g%ny = n2
        this%g%nn2min = 0
        this%g%nn2max = this%g%ny - 1

        this%g%n2min = this%g%nn2min
        this%g%n2max = this%g%nn2max
        this%nf2 = this%g%ny
        this%nf1min = this%g%nn1min
        this%nf2min = this%g%nn2min
        this%nf2max = this%g%nn2max
        this%g%m2min = this%g%nn2min
        this%g%m2max = this%g%nn2max
#endif
#if(_DIM_>=3)
        this%g%nz = n3
        this%g%nn3min = 0
        this%g%nn3max = this%g%nz - 1

        this%g%n3min = this%g%nn3min
        this%g%n3max = this%g%nn3max
        this%nf2 = this%g%nz
        this%nf3min = this%g%nn3min
        this%nf3max = this%g%nn3max
        this%g%m3min = this%g%nn3min
        this%g%m3max = this%g%nn3max
#endif  

#if(_DIM_==1)
        this%g%alloc_size = this%g%m1max-this%g%m1min+1
#elif(_DIM_==2)
        this%g%alloc_size = (this%g%m1max-this%g%m1min+1)*(this%g%m2max-this%g%m2min+1)
#elif(_DIM_==3)
        this%g%alloc_size = (this%g%m1max-this%g%m1min+1)*(this%g%m2max-this%g%m2min+1) &
                         *(this%g%m3max-this%g%m3min+1)
#endif

        allocate( this%g%nodes_x(this%g%n1min:this%g%n1max ) ) 
        allocate( this%g%weights_x(this%g%n1min:this%g%n1max ) ) 
        allocate( this%eigenvalues1(this%nf1min:this%nf1max ) ) 
        allocate( this%H1(this%g%n1min:this%g%n1max, this%nf1min:this%nf1max ) ) 
#if(_DIM_>=2)
        allocate( this%g%nodes_y(this%g%n2min:this%g%n2max ) ) 
        allocate( this%g%weights_y(this%g%n2min:this%g%n2max ) ) 
        allocate( this%eigenvalues2(this%nf2min:this%nf2max ) ) 
        allocate( this%H2(this%g%n2min:this%g%n2max, this%nf2min:this%nf2max ) ) 
#endif
#if(_DIM_>=3)
        allocate( this%g%nodes_z(this%g%n3min:this%g%n3max ) ) 
        allocate( this%g%weights_z(this%g%n3min:this%g%n3max ) ) 
        allocate( this%eigenvalues3(this%nf3min:this%nf3max ) ) 
        allocate( this%H3(this%g%n3min:this%g%n3max, this%nf3min:this%nf3max ) ) 
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
        class(S(tensorial)), intent(inout) :: this

        deallocate( this%g%nodes_x )
        deallocate( this%g%weights_x )
        deallocate( this%H1 )
        deallocate( this%eigenvalues1 )
#if(_DIM_>=2)
        deallocate( this%g%nodes_y )
        deallocate( this%g%weights_y )
        deallocate( this%H2 )
        deallocate( this%eigenvalues2 )
#endif
#if(_DIM_>=3)
        deallocate( this%g%nodes_z )
        deallocate( this%g%weights_z )
        deallocate( this%H3 )
        deallocate( this%eigenvalues3 )
#endif
#ifdef _OPENMP
        deallocate( this%g%jj )
        deallocate( this%jf )
#endif        
    end subroutine finalize_method



    function new_wf(m, u, coefficient) result(this)
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_ptr, c_loc
        type(S(wf_tensorial)) :: this
        class(S(tensorial)), target, intent(inout) :: m
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
    end function new_wf


    function clone(this) 
        class(S(wf_tensorial)), intent(inout) :: this
        class(_WAVE_FUNCTION_), pointer :: clone
        type(S(wf_tensorial)), pointer :: p

        allocate( p )
        p = S(wf_tensorial)(this%m, coefficient=this%coefficient)
        clone => p
    end function clone


    subroutine finalize_wf(this)
        use, intrinsic :: iso_c_binding, only: c_loc
        class(S(wf_tensorial)), intent(inout) :: this
        !call fftw_free(c_loc(this%up))
        call fftw_free(c_loc(this%up(1)))
    end subroutine finalize_wf


    subroutine set(this, f)
        class(S(wf_tensorial)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), external :: f
#ifdef _REAL_        
        call this%m%g%set_real_gridfun(this%u, f)
#else
        call this%m%g%set_complex_gridfun(this%u, f)
#endif        
        this%is_real_space = .true.
    end subroutine set


    subroutine set_t(this, f, t)
        class(S(wf_tensorial)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), external :: f
        real(kind=prec), intent(in) :: t
#ifdef _REAL_        
        call this%m%g%set_t_real_gridfun(this%u, f, t)
#else
        call this%m%g%set_t_complex_gridfun(this%u, f, t)
#endif        
        this%is_real_space = .true.
        this%time = t 
    end subroutine set_t


#ifndef _REAL_        
   subroutine rset(this, f)
        class(S(wf_tensorial)), intent(inout) :: this
        real(kind=prec), external :: f
        call this%m%g%rset_complex_gridfun(this%u, f)
        this%is_real_space = .true.
    end subroutine rset

   subroutine rset_t(this, f, t)
        class(S(wf_tensorial)), intent(inout) :: this
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t
        call this%m%g%rset_t_complex_gridfun(this%u, f, t)
        this%is_real_space = .true.
        this%time = t 
    end subroutine rset_t

#endif        

    function norm(this) 
        class(S(wf_tensorial)), intent(inout) :: this
        real(kind=prec) :: norm
        
        call this%to_real_space 
#ifdef _REAL_        
        norm = this%m%g%norm_real_gridfun(this%u)
#else
        norm = this%m%g%norm_complex_gridfun(this%u)
#endif        
    end function norm

    subroutine normalize(this, norm)
        class(S(wf_tensorial)), intent(inout) :: this
        real(kind=prec), intent(out), optional :: norm
        real(kind=prec) :: norm1
        norm1 = this%norm()
        if (present(norm)) then
            norm = norm1
        end if    
        this%u = (1.0_prec/norm1)*this%u
    end subroutine normalize

    subroutine scale(this, factor)
        class(S(wf_tensorial)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: factor
#ifdef _REAL_
        call this%m%g%scale_real_gridfun(this%u, factor)
#else        
        call this%m%g%scale_complex_gridfun(this%u, factor)
#endif        
    end subroutine scale

    subroutine axpy(this, other, factor)
         class(S(wf_tensorial)), intent(inout) :: this
         class(_WAVE_FUNCTION_), intent(inout) :: other
         _COMPLEX_OR_REAL_(kind=prec), intent(in) :: factor
        select type (other)
        class is (S(wf_tensorial))
        if (.not.associated(other%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if 
        call this%to_real_space 
        call other%to_real_space 
#ifdef _REAL_
        call this%m%g%axpy_real_gridfun(this%u, other%u, factor)
#else        
        call this%m%g%axpy_complex_gridfun(this%u, other%u, factor)
#endif        
        class default
           stop "E: wave functions not belonging to the same method"
        end select
    end subroutine axpy



    subroutine copy(this, source)
        class(S(wf_tensorial)), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: source 
        select type (source)
        class is (S(wf_tensorial))
        if (.not.associated(source%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    
       
        this%is_real_space = source%is_real_space
        this%time = source%time
        
#ifdef _REAL_
        call this%m%g%copy_real_gridfun(this%u, source%u)
#else        
        call this%m%g%copy_complex_gridfun(this%u, source%u)
#endif        
        class default
           stop "E: wave functions not belonging to the same method"
        end select
    end subroutine copy


    function distance(this, wf) result(d)
        class(S(wf_tensorial)), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: wf 
        real(kind=prec) :: d
        select type (wf)
        class is (S(wf_tensorial))
        if (.not.associated(wf%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if  
        if (associated(this%u, wf%u)) then
            d = 0.0_prec
            return
        end if

#ifdef _REAL_
        call this%axpy(wf, -1.0_prec)
        d = this%norm()
        call this%axpy(wf, +1.0_prec)
#else
        call this%axpy(wf, cmplx(-1.0_prec, 0.0_prec, kind=prec))
        d = this%norm()
        call this%axpy(wf, cmplx(-1.0_prec, 0.0_prec, kind=prec))
#endif        

        class default
           stop "E: wave functions not belonging to the same method"
        end select
    end function distance


    subroutine save(this, filename)
#ifdef _NO_HDF5_
        class(S(wf_tensorial)), intent(inout) :: this
        character(len=*), intent(in) :: filename
        print *, "W: save not implemented"
#else
#ifdef _QUADPRECISION_
        use tssmq_hdf5
#else
        use tssm_hdf5
#endif        
        class(S(wf_tensorial)), intent(inout) :: this
        character(len=*), intent(in) :: filename
        
        call this%to_real_space

#ifdef _REAL_        
        call hdf5_save_real_gridfun(this%m%g, this%u, filename, trim(this%m%dset_name))
#else        
        call hdf5_save_complex_gridfun(this%m%g, this%u, filename, &
                     trim(this%m%dset_name_real), trim(this%m%dset_name_imag))
#endif        

#if(_DIM_==1)        
        call hdf5_write_attributes(filename, (/ "grid" /), &
          (/ "tensorial_1D" /), &
          (/ "nx" /), (/ this%m%g%nx /), &
          (/ character(len=0) :: /), (/ real(prec) :: /) ) 
#elif(_DIM_==2)        
        call hdf5_write_attributes(filename, (/ "grid" /), &
          (/ "tensorial_2D" /), &
          (/ "nx", "ny" /), (/ this%m%g%nx, this%m%g%ny /), &
          (/ character(len=0) :: /), (/ real(prec) :: /) ) 
#elif(_DIM_==3)        
        call hdf5_write_attributes(filename, (/ "grid" /), &
          (/ "tensorial_3D" /), &
          (/ "nx", "ny", "nz" /), (/ this%m%g%nx, this%m%g%ny, this%m%g%nz /), &
          (/ character(len=0) :: /), (/ real(prec) :: /) ) 
#endif        
!TODO save nodes !!!
#endif        
    end subroutine save



    subroutine load(this, filename)
#ifdef _NO_HDF5_
        class(S(wf_tensorial)), intent(inout) :: this
        character(len=*), intent(in) :: filename
        print *, "W: load not implemented"
#else
#ifdef _QUADPRECISION_
        use tssmq_hdf5
#else
        use tssm_hdf5
#endif        
        class(S(wf_tensorial)), intent(inout) :: this
        character(len=*), intent(in) :: filename
#ifdef _QUADPRECISION_
        real(kind=prec), parameter :: eps = epsilon(1.0_8)
#else        
        real(kind=prec), parameter :: eps = epsilon(1.0_prec)
#endif
        character(len=20) :: s
        integer :: n, len


!TODO check compatibiliyt of grids!!!     
       call hdf5_read_string_attribute(filename, "grid", s, len)      
#if(_DIM_==1)
       if (s(1:len)/="tensorial_1D") then
              stop "E: wrong grid type"
       end if
#elif(_DIM_==2)
       if (s(1:len)/="tensorial_2D") then
              stop "E: wrong grid type"
       end if
#elif(_DIM_==3)
       if (s(1:len)/="tensorial_3D") then
              stop "E: wrong grid type"
       end if
#endif
       if (hdf5_read_integer_attribute(filename, "nx") /= this%m%g%nx) then
              stop "E: incompatible grids"
       end if       
#if(_DIM_>=2)
       if (hdf5_read_integer_attribute(filename, "ny") /= this%m%g%ny) then 
              stop "E: incompatible grids"
       end if       
#endif
#if(_DIM_>=3)
       if (hdf5_read_integer_attribute(filename, "nz") /= this%m%g%nz) then 
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


   subroutine to_real_space(this)
        class(S(wf_tensorial)), intent(inout) :: this
#if(_DIM_==3)        
        integer :: n1, n2, n3
#endif

        if (.not. this%is_real_space) then
!xxx$OMP PARALLEL WORKSHARE 
#if(_DIM_==1)        
           this%u = matmul(this%m%H1, this%uf)
#elif(_DIM_==2)        
           this%u = matmul(this%m%H1, this%uf)
           this%u = matmul(this%uf, transpose(this%m%H2))
#elif(_DIM_==3)        
           n1 = this%m%nf1max - this%m%nf1min + 1
           n2 = this%m%nf2max - this%m%nf2min + 1
           n3 = this%m%nf3max - this%m%nf3min + 1
           this%u = reshape( matmul(this%m%H1, reshape(this%uf, (/ n1, n2*n3 /) )), (/ n1, n2, n3 /) )
           this%u = reshape(reshape( matmul(this%m%H2,reshape( &
                            reshape(this%uf, (/ n2, n1, n3 /), order=(/ 2,1,3 /) ) , (/ n2, n1*n3 /) )), &
                            (/ n2, n1, n3 /) ), (/n1, n2, n3 /), order=(/ 2, 1, 3 /) )
           this%u = reshape( matmul(reshape(this%uf, (/ n1*n2, n3 /) ), transpose(this%m%H3)), (/ n1, n2, n3 /) ) 
#endif
!xxx$OMP END PARALLEL WORKSHARE 
           this%is_real_space = .true.
        end if     
    end subroutine to_real_space

    
    subroutine to_frequency_space(this)
        class(S(wf_tensorial)), intent(inout) :: this
#if(_DIM_==3)        
        integer :: n1, n2, n3
#endif

        if (this%is_real_space) then
!xxx$OMP PARALLEL WORKSHARE 
#if(_DIM_==1)        
           this%u = matmul(transpose(this%m%H1), this%m%g%weights_x*this%uf)
#elif(_DIM_==2)        
           this%u = matmul(transpose(this%m%H1),spread(this%m%g%weights_x, 2, this%m%g%n2max-this%m%g%n2min+1)*this%uf)
           this%u = matmul(spread(this%m%g%weights_y, 1, this%m%g%n1max-this%m%g%n1min+1)*this%uf, this%m%H2)
#elif(_DIM_==3)        
           n1 = this%m%nf1max - this%m%nf1min + 1
           n2 = this%m%nf2max - this%m%nf2min + 1
           n3 = this%m%nf3max - this%m%nf3min + 1
           this%u = reshape( matmul(transpose(this%m%H1), &
               reshape(this%uf, (/ n1, n2*n3 /) )*spread(this%m%g%weights_x, 2, n2*n3) &
               ), (/ n1, n2, n3 /) )
           this%u = reshape(reshape( matmul(transpose(this%m%H2), &
               reshape( reshape(this%uf, (/ n2, n1, n3 /), order=(/ 2,1,3 /) ) , (/ n2, n1*n3 /) ) &
                  * spread(this%m%g%weights_y, 2, n1*n3) &
               ),  (/ n2, n1, n3 /) ), (/ n1, n2, n3 /), order=(/ 2, 1, 3 /) ) 
           this%u = reshape( matmul( &
               reshape(this%uf, (/ n1*n2, n3 /) )*spread(this%m%g%weights_z, 1, n1*n2), &
               this%m%H3), (/ n1, n2, n3 /) ) 
#endif
!xxx$OMP END PARALLEL WORKSHARE 
           this%is_real_space = .false.
        end if     
    end subroutine to_frequency_space


    function norm_in_frequency_space(this) result(norm)
        class(S(wf_tensorial)), intent(inout) :: this
        real(kind=prec) :: norm

#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: norm1 
#endif
        
        call this%to_frequency_space 
!xxx$OMP PARALLEL WORKSHARE 
#ifdef _REAL_
        norm = sum( this%uf**2 )
#else
        norm = sum(real(this%uf,prec)**2 + aimag(this%uf)**2) 
#endif
!xxx$OMP END PARALLEL WORKSHARE 
#ifdef _MPI_
        call MPI_Reduce(norm, norm1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            norm = sqrt(norm1) 
        end if 
        call MPI_Bcast(norm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#else 
        norm = sqrt(norm)
#endif
    end function norm_in_frequency_space


    subroutine propagate_A(this, dt)
        class(S(wf_tensorial)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt
        if (dt==0.0_prec) then
            return
        end if  

        if (this%m%propagate_time_together_with_A) then
            call this%propagate_time(dt)
        end if

        call this%to_frequency_space

!xxx$OMP PARALLEL WORKSHARE 
#if(_DIM_==1)
        this%uf = exp((dt*this%coefficient)*this%m%eigenvalues1) * this%uf
#elif(_DIM_==2)
        this%uf = spread(exp((dt*this%coefficient)*this%m%eigenvalues1),2, this%m%nf2max-this%m%nf2min+1) * this%uf
        this%uf = spread(exp((dt*this%coefficient)*this%m%eigenvalues2),1, this%m%nf1max-this%m%nf1min+1) * this%uf
#elif(_DIM_==3)
        this%uf = spread(spread(exp((dt*this%coefficient)*this%m%eigenvalues1), 2, this%m%nf2max-this%m%nf2min+1), &
                                3, this%m%nf3max-this%m%nf3min+1) * this%uf
        this%uf = spread(spread(exp((dt*this%coefficient)*this%m%eigenvalues2), 1, this%m%nf1max-this%m%nf1min+1), &
                                3, this%m%nf3max-this%m%nf3min+1) * this%uf
        this%uf = spread(spread(exp((dt*this%coefficient)*this%m%eigenvalues3), 1, this%m%nf1max-this%m%nf1min+1), &
                                2, this%m%nf2max-this%m%nf2min+1) * this%uf
#endif
!xxx$OMP END PARALLEL WORKSHARE 
    end subroutine propagate_A


    subroutine propagate_A_derivative(this, wf, dt)
        class(S(wf_tensorial)), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: wf 
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt
        select type (wf)
        class is (S(wf_tensorial))
        if (.not.associated(wf%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    

        call wf%propagate_A(dt)
        call this%propagate_A(dt)

        class default
           stop "E: wrong spectral method for schroedinger wave function"
        end select
    end subroutine propagate_A_derivative


    subroutine add_apply_A(this, wf, coefficient)
        class(S(wf_tensorial)), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: wf 
        _COMPLEX_OR_REAL_(kind=prec), intent(in), optional :: coefficient
        _COMPLEX_OR_REAL_(kind=prec) :: C 

        select type (wf)
        class is (S(wf_tensorial))
        
        if (.not.associated(wf%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    

        C = this%coefficient
        if (present(coefficient)) then
            C = C*coefficient
        end if    

        call this%to_frequency_space
        call wf%to_frequency_space

!xxx$OMP PARALLEL WORKSHARE 
#if(_DIM_==1)
        wf%uf = wf%uf + C*this%m%eigenvalues1*this%uf
#elif(_DIM_==2)
        wf%uf = wf%uf + C*spread(this%m%eigenvalues1,2, this%m%nf2max-this%m%nf2min+1) * this%uf
        wf%uf = wf%uf + C*spread(this%m%eigenvalues2,1, this%m%nf1max-this%m%nf1min+1) * this%uf
#elif(_DIM_==3)
        wf%uf = wf%uf + C*spread(spread(this%m%eigenvalues1, 2, this%m%nf2max-this%m%nf2min+1), &
                                3, this%m%nf3max-this%m%nf3min+1) * this%uf
        wf%uf = wf%uf + C*spread(spread(this%m%eigenvalues2, 1, this%m%nf1max-this%m%nf1min+1), &
                                3, this%m%nf3max-this%m%nf3min+1) * this%uf
        wf%uf = wf%uf + C*spread(spread(this%m%eigenvalues3, 1, this%m%nf1max-this%m%nf1min+1), &
                                2, this%m%nf2max-this%m%nf2min+1) * this%uf
#endif
!xxx$OMP END PARALLEL WORKSHARE 
        class default
           stop "E: wave functions not belonging to the same method"
        end select
    end subroutine add_apply_A



#ifdef _QUADPRECISION_
end module S(tssmq_tensorial)
#else
end module S(tssm_tensorial)
#endif
