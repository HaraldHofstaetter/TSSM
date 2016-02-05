!!! TODO: OPENMP- Version of tssm_hermite_CPP.F90 !!!
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


module S(tssm_hermite)
    use tssm
    use tssm_grid
    use tssm_fourier_common   ! for fftw_alloc, etc.
    use tssm_hermite_common
    implicit none

    type, extends(spectral_method) :: S(hermite)
#if(_DIM_==1)
        type(grid_cartesian_1D) :: g 
#elif(_DIM_==2)
        type(grid_cartesian_2D) :: g 
#elif(_DIM_==3)
        type(grid_cartesian_3D) :: g 
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

        real(kind=prec) :: gamma_x 
        real(kind=prec), allocatable :: eigenvalues1(:)
        real(kind=prec), allocatable :: H_x(:,:)
#if(_DIM_>=2)
        real(kind=prec) :: gamma_y 
        real(kind=prec), allocatable :: eigenvalues2(:)
        real(kind=prec), allocatable :: H_y(:,:)
#endif 
#if(_DIM_>=3)
        real(kind=prec) :: gamma_z 
        real(kind=prec), allocatable :: eigenvalues3(:)
        real(kind=prec), allocatable :: H_z(:,:)
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
        procedure :: finalize => S(finalize_hermite)
        !final :: final_hermite_1D
        !! Fortran 2003 feature final seems to be not properly implemented
        !! in the gcc/gfortran compiler :(
    end type S(hermite)

    interface  S(hermite) ! constructor
        module procedure S(new_hermite)
    end interface S(hermite)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    type, extends(_WAVE_FUNCTION_) :: S(wf_hermite)
        class(S(hermite)), pointer :: m 

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
        procedure :: to_real_space => S(to_real_space_wf_hermite)
        procedure :: to_frequency_space => S(to_frequency_space_wf_hermite)
        procedure :: propagate_A => S(propagate_A_wf_hermite)
        procedure :: add_apply_A => S(add_apply_A_wf_hermite)
        procedure :: save => S(save_wf_hermite)
        procedure :: load => S(load_wf_hermite)
        procedure :: set => S(set_wf_hermite)
        procedure :: set_t => S(set_t_wf_hermite)
#ifndef _REAL_
        procedure :: rset => S(rset_wf_hermite)
        procedure :: rset_t => S(rset_t_wf_hermite)
#endif
        procedure :: clone => S(clone_wf_hermite)
        procedure :: finalize => S(finalize_wf_hermite)
        procedure :: copy => S(copy_wf_hermite)
        procedure :: norm2 => S(norm2_wf_hermite)
        procedure :: norm2_in_frequency_space => S(norm2_in_frequency_space_wf_hermite)
        procedure :: normalize => S(normalize_wf_hermite)
        procedure :: distance => S(distance_wf_hermite)
        procedure :: axpy => S(axpy_wf_hermite)
        procedure :: scale => S(scale_wf_hermite)
        !final :: S(final_wf_hermite)
    end type S(wf_hermite)

    interface S(wf_hermite) ! constructor
        module procedure S(new_wf_hermite)
    end interface S(wf_hermite)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains



#if(_DIM_==1)
    function S(new_hermite)(nx, gamma_x ) result(this)
#elif(_DIM_==2)
    function S(new_hermite)(nx, gamma_x, ny, gamma_y) result(this)
#elif(_DIM_==3)
    function S(new_hermite)(nx, gamma_x, ny, gamma_y, nz, gamma_z) result(this)
#endif
        type(S(hermite)) :: this
        integer, intent(in) :: nx
        real(kind=prec), intent(in) :: gamma_x 
#if(_DIM_>=2)
        integer, intent(in) :: ny
        real(kind=prec), intent(in) :: gamma_y 
#endif
#if(_DIM_>=3)
        integer, intent(in) :: nz
        real(kind=prec), intent(in) :: gamma_z 
#endif

        integer :: j

!adjust grid parameters, TODO MPI!!!


        this%g%nx = nx
        this%g%nn1min = 0
        this%g%nn1max = this%g%nx

        this%g%n1min = this%g%nn1min
        this%g%n1max = this%g%nn1max
        this%nf1min = this%g%nn1min
        this%nf1max = this%g%nn1max
        this%g%m1min = this%g%nn1min
        this%g%m1max = this%g%nn1max
#if(_DIM_>=2)
        this%g%ny = ny
        this%g%nn2min = 0
        this%g%nn2max = this%g%ny

        this%g%n2min = this%g%nn2min
        this%g%n2max = this%g%nn2max
        this%nf2min = this%g%nn2min
        this%nf2max = this%g%nn2max
        this%g%m2min = this%g%nn2min
        this%g%m2max = this%g%nn2max
#endif
#if(_DIM_>=3)
        this%g%nz = nz
        this%g%nn3min = 0
        this%g%nn3max = this%g%nz

        this%g%n3min = this%g%nn3min
        this%g%n3max = this%g%nn3max
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

        this%gamma_x = gamma_x
        allocate( this%g%nodes_x(this%g%n1min:this%g%n1max ) ) 
        allocate( this%g%weights_x(this%g%n1min:this%g%n1max ) ) 
        allocate( this%eigenvalues1(this%nf1min:this%nf1max ) ) 
        allocate( this%H_x(this%g%n1min:this%g%n1max, this%nf1min:this%nf1max ) ) 
        call hermite_scaled_coeffs(nx, this%g%nodes_x,  this%g%weights_x, this%H_x, this%gamma_x)
        this%eigenvalues1 = (/ ( -this%gamma_x*(0.5_prec + real(j, prec)), j = this%nf1min, this%nf1max) /)
#if(_DIM_>=2)
        this%gamma_y = gamma_y
        allocate( this%g%nodes_y(this%g%n2min:this%g%n2max ) ) 
        allocate( this%g%weights_y(this%g%n2min:this%g%n2max ) ) 
        allocate( this%eigenvalues2(this%nf2min:this%nf2max ) ) 
        allocate( this%H_y(this%g%n2min:this%g%n2max, this%nf2min:this%nf2max ) ) 
        call hermite_scaled_coeffs(ny, this%g%nodes_y,  this%g%weights_y, this%H_y, this%gamma_y)
        this%eigenvalues2 = (/ ( -this%gamma_y*(0.5_prec + real(j, prec)), j = this%nf2min, this%nf2max) /)
#endif
#if(_DIM_>=3)
        this%gamma_z = gamma_z
        allocate( this%g%nodes_z(this%g%n3min:this%g%n3max ) ) 
        allocate( this%g%weights_z(this%g%n3min:this%g%n3max ) ) 
        allocate( this%eigenvalues3(this%nf3min:this%nf3max ) ) 
        allocate( this%H_z(this%g%n3min:this%g%n3max, this%nf3min:this%nf3max ) ) 
        call hermite_scaled_coeffs(nz, this%g%nodes_z,  this%g%weights_z, this%H_z, this%gamma_z)
        this%eigenvalues3 = (/ ( -this%gamma_z*(0.5_prec + real(j, prec)), j = this%nf3min, this%nf3max) /)
#endif

    end function S(new_hermite)

    subroutine S(finalize_hermite)(this)
        class(S(hermite)), intent(inout) :: this

        deallocate( this%g%nodes_x )
        deallocate( this%g%weights_x )
        deallocate( this%H_x )
        deallocate( this%eigenvalues1 )
#if(_DIM_>=2)
        deallocate( this%g%nodes_y )
        deallocate( this%g%weights_y )
        deallocate( this%H_y )
        deallocate( this%eigenvalues2 )
#endif
#if(_DIM_>=3)
        deallocate( this%g%nodes_z )
        deallocate( this%g%weights_z )
        deallocate( this%H_z )
        deallocate( this%eigenvalues3 )
#endif
#ifdef _OPENMP
        deallocate( this%g%jj )
        deallocate( this%jf )
#endif        
    end subroutine S(finalize_hermite)



    function S(new_wf_hermite)(m, u, coefficient) result(this)
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_ptr, c_loc
        type(S(wf_hermite)) :: this
        class(S(hermite)), target, intent(inout) :: m
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
    end function S(new_wf_hermite)


    function S(clone_wf_hermite)(this) result(clone)
        class(S(wf_hermite)), intent(inout) :: this
        class(_WAVE_FUNCTION_), pointer :: clone
        type(S(wf_hermite)), pointer :: p

        allocate( p )
        p = S(new_wf_hermite)(this%m, coefficient=this%coefficient)
        clone => p
    end function S(clone_wf_hermite)


    subroutine S(finalize_wf_hermite)(this)
        use, intrinsic :: iso_c_binding, only: c_loc
        class(S(wf_hermite)), intent(inout) :: this
        !call fftw_free(c_loc(this%up))
        call fftw_free(c_loc(this%up(1)))
    end subroutine S(finalize_wf_hermite)



    subroutine S(set_wf_hermite)(this, f)
        class(S(wf_hermite)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), external :: f
#ifdef _REAL_        
        call this%m%g%set_real_gridfun(this%u, f)
#else
        call this%m%g%set_complex_gridfun(this%u, f)
#endif        
        this%is_real_space = .true.
    end subroutine S(set_wf_hermite)

    subroutine S(set_t_wf_hermite)(this, f, t)
        class(S(wf_hermite)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), external :: f
        real(kind=prec), intent(in) :: t
#ifdef _REAL_        
        call this%m%g%set_t_real_gridfun(this%u, f, t)
#else
        call this%m%g%set_t_complex_gridfun(this%u, f, t)
#endif        
        this%is_real_space = .true.
    end subroutine S(set_t_wf_hermite)


#ifndef _REAL_        
   subroutine S(rset_wf_hermite)(this, f)
        class(S(wf_hermite)), intent(inout) :: this
        real(kind=prec), external :: f
        call this%m%g%rset_complex_gridfun(this%u, f)
        this%is_real_space = .true.
    end subroutine S(rset_wf_hermite)

   subroutine S(rset_t_wf_hermite)(this, f, t)
        class(S(wf_hermite)), intent(inout) :: this
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t
        call this%m%g%rset_t_complex_gridfun(this%u, f, t)
        this%is_real_space = .true.
    end subroutine S(rset_t_wf_hermite)

#endif        

    function S(norm2_wf_hermite)(this) result(n)
        class(S(wf_hermite)), intent(inout) :: this
        real(kind=prec) :: n
        
        call this%to_real_space 
#ifdef _REAL_        
        n = this%m%g%norm2_real_gridfun(this%u)
#else
        n = this%m%g%norm2_complex_gridfun(this%u)
#endif        
    end function S(norm2_wf_hermite)

    subroutine S(normalize_wf_hermite)(this, norm)
        class(S(wf_hermite)), intent(inout) :: this
        real(kind=prec), intent(out), optional :: norm
        real(kind=prec) :: n
        n = this%norm2()
        if (present(norm)) then
            norm = n
        end if    
        this%u = (1.0_prec/n)*this%u
    end subroutine S(normalize_wf_hermite)

    subroutine S(scale_wf_hermite)(this, factor)
        class(S(wf_hermite)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: factor
#ifndef _REAL_
        if(aimag(factor)==0.0_prec) then
            this%u = real(factor, kind=prec)*this%u
        else
#endif        
            this%u = factor*this%u
#ifndef _REAL_
        end if
#endif        
    end subroutine S(scale_wf_hermite)

    subroutine S(axpy_wf_hermite)(this, other, factor)
         class(S(wf_hermite)), intent(inout) :: this
         class(_WAVE_FUNCTION_), intent(inout) :: other
         _COMPLEX_OR_REAL_(kind=prec), intent(in) :: factor
        select type (other)
        class is (S(wf_hermite))
        if (.not.associated(other%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if 
        call this%to_real_space 
        call other%to_real_space 
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
        class default
           stop "E: wave functions not belonging to the same method"
        end select
    end subroutine S(axpy_wf_hermite)



    subroutine S(copy_wf_hermite)(this, source)
        class(S(wf_hermite)), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: source 
        select type (source)
        class is (S(wf_hermite))
        if (.not.associated(source%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    
       
        this%is_real_space = source%is_real_space

        if (this%is_real_space) then
            this%u = source%u
        else    
            this%uf = source%uf 
        end if

        class default
           stop "E: wave functions not belonging to the same method"
        end select
    end subroutine S(copy_wf_hermite)


    function S(distance_wf_hermite)(this, wf) result(n)
        class(S(wf_hermite)), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: wf 
        real(kind=prec) :: n
        select type (wf)
        class is (S(wf_hermite))
        if (.not.associated(wf%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if  
        if (associated(this%u, wf%u)) then
            n = 0.0_prec
            return
        end if
  
        !!! TODO handle norm in frequency space without transforming
#if 1
        call this%to_real_space 
        call wf%to_real_space 
#ifdef _REAL_        
        n = this%m%g%norm2_real_gridfun(this%u-wf%u)
#else
        n = this%m%g%norm2_complex_gridfun(this%u-wf%u)
#endif      
#else
        call this%to_frequency_space 
        call wf%to_frequency_space 
!$OMP PARALLEL WORKSHARE 
        this%uf = this%uf - wf%uf
!$OMP END PARALLEL WORKSHARE 
        n = this%norm2()
!$OMP PARALLEL WORKSHARE 
        this%uf = this%uf + wf%uf
!$OMP END PARALLEL WORKSHARE 
#endif
        class default
           stop "E: wave functions not belonging to the same method"
        end select
    end function S(distance_wf_hermite)



    subroutine S(save_wf_hermite)(this, filename)
#ifdef _NO_HDF5_
        class(S(wf_hermite)), intent(inout) :: this
        character(len=*), intent(in) :: filename
        print *, "W: save not implemented"
#else
        use tssm_hdf5
        class(S(wf_hermite)), intent(inout) :: this
        character(len=*), intent(in) :: filename
        
        call this%to_real_space

#ifdef _REAL_        
        call hdf5_save_real_gridfun(this%m%g, this%u, filename, trim(this%m%dset_name))
#else        
        call hdf5_save_complex_gridfun(this%m%g, this%u, filename, &
                     trim(this%m%dset_name_real), trim(this%m%dset_name_imag))
#endif        

#if(_DIM_==1)        
        call hdf5_write_attributes(filename, (/ "grid      ", "nodes_type" /), &
          (/ "cartesian_1D", "hermite     " /), &
          (/ "nx" /), (/ this%m%g%nx /), &
          (/ "gamma_x" /), (/ this%m%gamma_x /) ) 
#elif(_DIM_==2)        
        call hdf5_write_attributes(filename, (/ "grid      ", "nodes_type" /), &
          (/ "cartesian_2D", "hermite     " /), &
          (/ "nx", "ny" /), (/ this%m%g%nx, this%m%g%ny /), &
          (/ "gamma_x", "gamma_y" /), (/ this%m%gamma_x,  this%m%gamma_y /) ) 
#elif(_DIM_==3)        
        call hdf5_write_attributes(filename, (/ "grid      ", "nodes_type" /), &
          (/ "cartesian_3D", "hermite     " /), &
          (/ "nx", "ny", "nz" /), (/ this%m%g%nx, this%m%g%ny, this%m%g%nz /), &
          (/ "gamma_x", "gamma_y", "gamma_z" /), (/ this%m%gamma_x, this%m%gamma_y, this%m%gamma_z /) ) 
#endif        
!TODO save nodes !!!
#endif        
    end subroutine S(save_wf_hermite)



    subroutine S(load_wf_hermite)(this, filename)
#ifdef _NO_HDF5_
        class(S(wf_hermite)), intent(inout) :: this
        character(len=*), intent(in) :: filename
        print *, "W: load not implemented"
#else
        use tssm_hdf5
        class(S(wf_hermite)), intent(inout) :: this
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
       if (s(1:len)/="cartesian_1D") then
              stop "E: wrong grid type"
       end if
#elif(_DIM_==2)
       if (s(1:len)/="cartesian_2D") then
              stop "E: wrong grid type"
       end if
#elif(_DIM_==3)
       if (s(1:len)/="cartesian_3D") then
              stop "E: wrong grid type"
       end if
#endif
       call hdf5_read_string_attribute(filename, "nodes_type", s, len)      
       if (s(1:len)/="hermite") then
              stop "E: wrong nodes type"
       end if
       
       if ((hdf5_read_integer_attribute(filename, "nx") /= this%m%g%nx) .or. &
           (abs(hdf5_read_double_attribute(filename, "gamma_x")- this%m%gamma_x)>eps)) then
              stop "E: incompatible grids"
       end if       
#if(_DIM_>=2)
       if ((hdf5_read_integer_attribute(filename, "ny") /= this%m%g%ny) .or. &
           (abs(hdf5_read_double_attribute(filename, "gamma_y")- this%m%gamma_y)>eps)) then
              stop "E: incompatible grids"
       end if       
#endif
#if(_DIM_>=3)
       if ((hdf5_read_integer_attribute(filename, "nz") /= this%m%g%nz) .or. &
           (abs(hdf5_read_double_attribute(filename, "gamma_z")- this%m%gamma_z)>eps)) then
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
    end subroutine S(load_wf_hermite)


   subroutine S(to_real_space_wf_hermite)(this)
        class(S(wf_hermite)), intent(inout) :: this
#if(_DIM_==3)        
        integer :: nx, ny, nz
#endif

        if (.not. this%is_real_space) then
!$OMP PARALLEL WORKSHARE 
#if(_DIM_==1)        
           this%u = matmul(this%m%H_x, this%uf)
#elif(_DIM_==2)        
           this%u = matmul(this%m%H_x, this%uf)
           this%u = matmul(this%uf, transpose(this%m%H_y))
#elif(_DIM_==3)        
           nx = this%m%nf1max - this%m%nf1min + 1
           ny = this%m%nf2max - this%m%nf2min + 1
           nz = this%m%nf3max - this%m%nf3min + 1
           this%u = reshape( matmul(this%m%H_x, reshape(this%uf, (/ nx, ny*nz /) )), (/ nx, ny, nz /) )
           this%u = reshape(reshape( matmul(this%m%H_y,reshape( &
                            reshape(this%uf, (/ ny, nx, nz /), order=(/ 2,1,3 /) ) , (/ ny, nx*nz /) )), &
                            (/ ny, nx, nz /) ), (/nx, ny, nz /), order=(/ 2, 1, 3 /) )
           this%u = reshape( matmul(reshape(this%uf, (/ nx*ny, nz /) ), transpose(this%m%H_z)), (/ nx, ny, nz /) ) 
#endif
!$OMP END PARALLEL WORKSHARE 
           this%is_real_space = .true.
        end if     
    end subroutine S(to_real_space_wf_hermite)

    
    subroutine S(to_frequency_space_wf_hermite)(this)
        class(S(wf_hermite)), intent(inout) :: this
#if(_DIM_==3)        
        integer :: nx, ny, nz
#endif

        if (this%is_real_space) then
!$OMP PARALLEL WORKSHARE 
#if(_DIM_==1)        
           this%u = matmul(transpose(this%m%H_x), this%m%g%weights_x*this%uf)
#elif(_DIM_==2)        
           this%u = matmul(transpose(this%m%H_x),spread(this%m%g%weights_x, 2, this%m%g%n2max-this%m%g%n2min+1)*this%uf)
           this%u = matmul(spread(this%m%g%weights_y, 1, this%m%g%n1max-this%m%g%n1min+1)*this%uf, this%m%H_y)
#elif(_DIM_==3)        
           nx = this%m%nf1max - this%m%nf1min + 1
           ny = this%m%nf2max - this%m%nf2min + 1
           nz = this%m%nf3max - this%m%nf3min + 1
           this%u = reshape( matmul(transpose(this%m%H_x), &
               reshape(this%uf, (/ nx, ny*nz /) )*spread(this%m%g%weights_x, 2, ny*nz) &
               ), (/ nx, ny, nz /) )
           this%u = reshape(reshape( matmul(transpose(this%m%H_y), &
               reshape( reshape(this%uf, (/ ny, nx, nz /), order=(/ 2,1,3 /) ) , (/ ny, nx*nz /) ) &
                  * spread(this%m%g%weights_y, 2, nx*nz) &
               ),  (/ ny, nx, nz /) ), (/ nx, ny, nz /), order=(/ 2, 1, 3 /) ) 
           this%u = reshape( matmul( &
               reshape(this%uf, (/ nx*ny, nz /) )*spread(this%m%g%weights_z, 1, nx*ny), &
               this%m%H_z), (/ nx, ny, nz /) ) 
#endif
!$OMP END PARALLEL WORKSHARE 
           this%is_real_space = .false.
        end if     
    end subroutine S(to_frequency_space_wf_hermite)


    function S(norm2_in_frequency_space_wf_hermite)(this) result(N)
        class(S(wf_hermite)), intent(inout) :: this
        real(kind=prec) :: N

#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: n1 
#endif
        
        call this%to_frequency_space 
!$OMP PARALLEL WORKSHARE 
#ifdef _REAL_
        n = sum( this%uf**2 )
#else
        n = sum(real(this%uf,prec)**2 + aimag(this%uf)**2) 
#endif
!$OMP END PARALLEL WORKSHARE 
#ifdef _MPI_
        call MPI_Reduce(n, n1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            n = sqrt(n1) 
        end if 
        call MPI_Bcast(N, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#else 
        n = sqrt(n)
#endif
    end function S(norm2_in_frequency_space_wf_hermite)


    subroutine S(propagate_A_wf_hermite)(this, dt)
        class(S(wf_hermite)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt
        call this%to_frequency_space

!$OMP PARALLEL WORKSHARE 
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
!$OMP END PARALLEL WORKSHARE 
    end subroutine S(propagate_A_wf_hermite)



    subroutine S(add_apply_A_wf_hermite)(this, wf, coefficient)
        class(S(wf_hermite)), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: wf 
        _COMPLEX_OR_REAL_(kind=prec), intent(in), optional :: coefficient
        _COMPLEX_OR_REAL_(kind=prec) :: C 

        select type (wf)
        class is (S(wf_hermite))
        
        if (.not.associated(wf%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    

        C = this%coefficient
        if (present(coefficient)) then
            C = C*coefficient
        end if    

        call this%to_frequency_space
        call wf%to_frequency_space

!$OMP PARALLEL WORKSHARE 
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
!$OMP END PARALLEL WORKSHARE 
        class default
           stop "E: wave functions not belonging to the same method"
        end select
    end subroutine S(add_apply_A_wf_hermite)



end module S(tssm_hermite)
