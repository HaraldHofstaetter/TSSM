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
module S(tssmq_multicomponent_fourier)    
    use S(tssmq_fourier)
#else
module S(tssm_multicomponent_fourier)    
    use S(tssm_fourier)
#endif    
    implicit none

    type, extends(spectral_method) :: S(multicomponent_fourier)
        class(S(fourier)), pointer :: b => null() ! type of components
        integer :: nc                             ! number of components
        _COMPLEX_OR_REAL_(kind=prec), allocatable :: coefficients(:)
#ifdef _REAL_
       character(len=32), allocatable :: dset_names(:)
#else
       character(len=32), allocatable :: dset_names_real(:)
       character(len=32), allocatable :: dset_names_imag(:)
#endif
        integer :: method_for_B = 2 
    end type S(multicomponent_fourier)

    interface S(multicomponent_fourier)
        module procedure :: S(new_multicomponent_fourier)
    end interface S(multicomponent_fourier)

    type :: S(pointer_to_wf_fourier)
        class(S(wf_fourier)), pointer :: p
    end type S(pointer_to_wf_fourier)


    type, extends(_WAVE_FUNCTION_) :: S(wf_multicomponent_fourier)
        class(S(multicomponent_fourier)), pointer :: m 
        class(S(pointer_to_wf_fourier)), allocatable :: C(:)  ! .. components
    contains 
        procedure :: to_real_space => S(to_real_space_wf_multicomponent_fourier)
        procedure :: to_frequency_space => S(to_frequency_space_wf_multicomponent_fourier)
        procedure :: propagate_A => S(propagate_A_wf_multicomponent_fourier)
        procedure :: load => S(load_wf_multicomponent_fourier)
        procedure :: save => S(save_wf_multicomponent_fourier)
        procedure :: set => S(set_wf_multicomponent_fourier)
        procedure :: set_t => S(set_t_wf_multicomponent_fourier)
#ifndef _REAL_
        procedure :: rset => S(rset_wf_multicomponent_fourier)
        procedure :: rset_t => S(rset_t_wf_multicomponent_fourier)
#endif
        procedure :: clone => S(clone_wf_multicomponent_fourier)
        procedure :: finalize => S(finalize_wf_multicomponent_fourier)
        procedure :: copy => S(copy_wf_multicomponent_fourier)
        procedure :: norm => S(norm_wf_multicomponent_fourier)
        procedure :: normalize => S(normalize_wf_multicomponent_fourier)
        procedure :: distance => S(distance_wf_multicomponent_fourier)
        procedure :: axpy => S(axpy_wf_multicomponent_fourier)
        procedure :: scale => S(scale_wf_multicomponent_fourier)
    end type S(wf_multicomponent_fourier)

    interface S(wf_multicomponent_fourier)
        module procedure S(new_wf_multicomponent_fourier)
    end interface S(wf_multicomponent_fourier)

contains

#if(_DIM_==1)                           
    function S(new_multicomponent_fourier)(coefficients, nx, xmin, xmax, transformation_type) result(this)
#elif(_DIM_==2)
    function S(new_multicomponent_fourier)(coefficients, nx, xmin, xmax, ny, ymin, ymax, &
                               transformation_type) result(this)
#elif(_DIM_==3)
    function S(new_multicomponent_fourier)(coefficients, nx, xmin, xmax, ny, ymin, ymax, &
                               nz, zmin, zmax, transformation_type) result(this)
#endif
        type(S(multicomponent_fourier)) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: coefficients(:)
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
        integer, intent(in), optional :: transformation_type

        character(len=4) :: s
        integer :: k
        type(S(fourier)), pointer :: b

        this%nc = ubound(coefficients, 1) 
        allocate( this%coefficients( this%nc ) )
        this%coefficients = coefficients
#ifdef _REAL_
        allocate( this%dset_names( this%nc ) )
        ! default names
        do k = 1, this%nc
            write (s,'(I2.2)') k
            this%dset_names(k) = "psi" // s 
        end do
#else        
        allocate( this%dset_names_real( this%nc ) )
        allocate( this%dset_names_imag( this%nc ) )
        ! default names
        do k = 1, this%nc
            write (s,'(I2.2)') k
            this%dset_names_real(k) = "psi" // s // "_real"
            this%dset_names_imag(k) = "psi" // s // "_imag"
        end do
#endif        
        if(.not.associated(this%b)) then
            allocate( b )
#if(_DIM_==1)
            b = S(fourier)(nx, xmin, xmax, transformation_type)
#elif(_DIM_==2)
            b = S(fourier)(nx, xmin, xmax, ny, ymin, ymax, transformation_type)
#elif(_DIM_==3)
            b = S(fourier)(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax, transformation_type)
#endif
            this%b => b
        end if
    end function S(new_multicomponent_fourier)


    function S(new_wf_multicomponent_fourier)(m) result(this)
        type(S(wf_multicomponent_fourier)) :: this
!        class(S(wf_multicomponent_fourier)) :: this
        class(S(multicomponent_fourier)), target, intent(inout) :: m
        integer :: j

        type(S(wf_fourier)), pointer :: p

        this%m => m
        allocate( this%C( m%nc ) )
        select type (b=>m%b)
        class is (S(fourier))
        do j = 1, m%nc
            !this%C(j) = S(new_wf_fourier)(m%S(fourier), coefficient=m%coefficients(j))
            allocate(  p )
            p = S(new_wf_fourier)(b, coefficient=m%coefficients(j))
            this%C(j)%p => p
        end do    
        end select
    end function S(new_wf_multicomponent_fourier)


   function S(clone_wf_multicomponent_fourier)(this) result(clone)
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        class(_WAVE_FUNCTION_), pointer :: clone
        type(S(wf_multicomponent_fourier)), pointer :: p
        
        allocate( p )
        p = S(wf_multicomponent_fourier)(this%m)
        clone => p
    end function S(clone_wf_multicomponent_fourier)


    subroutine S(finalize_wf_multicomponent_fourier)(this)
        use, intrinsic :: iso_c_binding, only: c_loc
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        integer :: j

        do j = 1, this%m%nc 
            call this%C(j)%p%finalize()
            deallocate( this%C(j)%p ) 
        end do   
        deallocate( this%C ) 
    end subroutine S(finalize_wf_multicomponent_fourier)

    subroutine S(copy_wf_multicomponent_fourier)(this, source)
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: source 
        integer :: j

        select type (source)
        class is (S(wf_multicomponent_fourier))
        if (.not.associated(source%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if   
        do j = 1, this%m%nc 
            call this%C(j)%p%copy(source%C(j)%p)
        end do   
        class default
           stop "E: wave functions not belonging to the same method"
        end select
    end subroutine S(copy_wf_multicomponent_fourier)



    subroutine S(to_real_space_wf_multicomponent_fourier)(this)
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        integer :: j
        do j = 1, this%m%nc 
            call this%C(j)%p%to_real_space()
        end do    
    end subroutine S(to_real_space_wf_multicomponent_fourier)

    subroutine S(to_frequency_space_wf_multicomponent_fourier)(this)
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        integer :: j
        do j = 1, this%m%nc 
            call this%C(j)%p%to_frequency_space()
        end do    
    end subroutine S(to_frequency_space_wf_multicomponent_fourier)




    subroutine S(propagate_A_wf_multicomponent_fourier)(this, dt)
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt
        integer :: j
        do j = 1, this%m%nc 
            call this%C(j)%p%propagate_A(dt)
        end do    
    end subroutine S(propagate_A_wf_multicomponent_fourier)


    function S(norm_wf_multicomponent_fourier)(this) result(n)
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        real(kind=prec) :: n
        integer :: j
        n = 0.0_prec 
        do j=1, this%m%nc
             n = n + ( this%C(j)%p%norm() )**2 
        end do
        n = sqrt(n)
    end function S(norm_wf_multicomponent_fourier)


    function S(distance_wf_multicomponent_fourier)(this, wf) result(n)
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: wf 
        real(kind=prec) :: n
        integer :: j
 
        select type (wf)
        class is (S(wf_multicomponent_fourier))
        
        if (.not.associated(wf%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    

        n = 0.0_prec 
        do j=1, this%m%nc 
             n = n + ( this%C(j)%p%distance(wf%C(j)%p) )**2 
        end do
        n = sqrt(n)

        class default
           stop "E: wave functions not belonging to the same method"
        end select
    end function S(distance_wf_multicomponent_fourier)


    subroutine S(scale_wf_multicomponent_fourier)(this, factor)
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: factor

        integer :: j

        do j = 1, this%m%nc 
            call this%C(j)%p%scale(factor)
        end do    
    end subroutine S(scale_wf_multicomponent_fourier)

    subroutine S(normalize_wf_multicomponent_fourier)(this, norm)
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        real(kind=prec), intent(out), optional :: norm
        real(kind=prec) :: n

        n = this%norm()
#ifdef _REAL_        
        call this%scale(1.0_prec/n)
#else
        call this%scale(cmplx(1.0_prec/n, kind=prec))
#endif
        if(present(norm)) then
           norm = n
        end if   
    end subroutine S(normalize_wf_multicomponent_fourier)


    subroutine S(axpy_wf_multicomponent_fourier)(this, other, factor)
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: other
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: factor

        integer :: j

        select type (other)
        class is (S(wf_multicomponent_fourier))
        if (.not.associated(other%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    

        do j = 1, this%m%nc 
            call this%C(j)%p%axpy(other%C(j)%p, factor)
        end do    

        class default
           stop "E: wave functions not belonging to the same method"
        end select
    end subroutine S(axpy_wf_multicomponent_fourier)



    subroutine S(set_wf_multicomponent_fourier)(this, f)
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), external :: f
        integer :: j
        do j=1, this%m%nc
            call this%C(j)%p%set(f_component)
        end do    
    contains

#if(_DIM_==1)        
        function f_component(x)
            _COMPLEX_OR_REAL_(kind=prec) :: f_component
            _COMPLEX_OR_REAL_(kind=prec), intent(in) :: x
            f_component = f(x, j)
        end function f_component
#elif(_DIM_==2)        
        function f_component(x, y)
            _COMPLEX_OR_REAL_(kind=prec) :: f_component
            _COMPLEX_OR_REAL_(kind=prec), intent(in) :: x, y
            f_component = f(x, y, j)
        end function f_component
#elif(_DIM_==3)        
        function f_component(x, y, z)
            _COMPLEX_OR_REAL_(kind=prec) :: f_component
            _COMPLEX_OR_REAL_(kind=prec), intent(in) :: x, y, z
            f_component = f(x, y, z, j)
        end function f_component
#endif        
    end subroutine S(set_wf_multicomponent_fourier)

    subroutine S(set_t_wf_multicomponent_fourier)(this, f, t)
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), external :: f
        real(kind=prec), intent(in) :: t
        integer :: j
        do j=1, this%m%nc
            call this%C(j)%p%set_t(f_component, t)
        end do    
    contains

#if(_DIM_==1)        
        function f_component(x, t)
            _COMPLEX_OR_REAL_(kind=prec) :: f_component
            _COMPLEX_OR_REAL_(kind=prec), intent(in) :: x
            real(kind=prec), intent(in) :: t
            f_component = f(x, t, j)
        end function f_component
#elif(_DIM_==2)        
        function f_component(x, y, t)
            _COMPLEX_OR_REAL_(kind=prec) :: f_component
            _COMPLEX_OR_REAL_(kind=prec), intent(in) :: x, y
            real(kind=prec), intent(in) :: t
            f_component = f(x, y, t, j)
        end function f_component
#elif(_DIM_==3)        
        function f_component(x, y, z, t)
            _COMPLEX_OR_REAL_(kind=prec) :: f_component
            _COMPLEX_OR_REAL_(kind=prec), intent(in) :: x, y, z
            real(kind=prec), intent(in) :: t
            f_component = f(x, y, z, t, j)
        end function f_component
#endif        
    end subroutine S(set_t_wf_multicomponent_fourier)




#ifndef _REAL_        
   subroutine S(rset_wf_multicomponent_fourier)(this, f)
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        real(kind=prec), external :: f
        integer :: j
        do j=1, this%m%nc
            call this%C(j)%p%rset(f_component)
        end do    
    contains

#if(_DIM_==1)        
        function f_component(x)
            real(kind=prec) :: f_component
            real(kind=prec), intent(in) :: x
            f_component = f(x, j)
        end function f_component
#elif(_DIM_==2)        
        function f_component(x, y)
            real(kind=prec) :: f_component
            real(kind=prec), intent(in) :: x, y
            f_component = f(x, y, j)
        end function f_component
#elif(_DIM_==3)        
        function f_component(x, y, z)
            real(kind=prec) :: f_component
            real(kind=prec), intent(in) :: x, y, z
            f_component = f(x, y, z, j)
        end function f_component
#endif        
    end subroutine S(rset_wf_multicomponent_fourier)


   subroutine S(rset_t_wf_multicomponent_fourier)(this, f, t)
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t
        integer :: j
        do j=1, this%m%nc
            call this%C(j)%p%rset_t(f_component, t)
        end do    
    contains

#if(_DIM_==1)        
        function f_component(x, t)
            real(kind=prec) :: f_component
            real(kind=prec), intent(in) :: x
            real(kind=prec), intent(in) :: t
            f_component = f(x, t, j)
        end function f_component
#elif(_DIM_==2)        
        function f_component(x, y, t)
            real(kind=prec) :: f_component
            real(kind=prec), intent(in) :: x, y
            real(kind=prec), intent(in) :: t
            f_component = f(x, y, t, j)
        end function f_component
#elif(_DIM_==3)        
        function f_component(x, y, z, t)
            real(kind=prec) :: f_component
            real(kind=prec), intent(in) :: x, y, z
            real(kind=prec), intent(in) :: t
            f_component = f(x, y, z, t, j)
        end function f_component
#endif        
    end subroutine S(rset_t_wf_multicomponent_fourier)

#endif        

    subroutine S(load_wf_multicomponent_fourier)(this, filename)
#ifdef _NO_HDF5_
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        character(len=*), intent(in) :: filename
        print *, "W: load not implemented"
#else
#ifdef _QUADPRECISION_
        use tssmq_hdf5
#else
        use tssm_hdf5
#endif        
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer :: k
#ifdef _QUADPRECISION_        
        real(kind=prec), parameter :: eps = epsilon(1.0_8)
#else        
        real(kind=prec), parameter :: eps = epsilon(1.0_prec)
#endif        
#if(_DIM_==1)
        type(grid_equidistant_1D) :: g 
#elif(_DIM_==2)
        type(grid_equidistant_2D) :: g 
#elif(_DIM_==3)
        type(grid_equidistant_3D) :: g 
#endif 
#ifdef _REAL_
        call hdf5_read_grid_attributes(g, filename, trim(this%m%b%dset_name))
#else        
        call hdf5_read_grid_attributes(g, filename, trim(this%m%b%dset_name_real))
#endif        
        !TODO do not duplicate this code segment from tssm_fourier_CPP.F90
        !TODO handle certain kinds of 'incompatible' grids
#if(_DIM_==1)
        if (.not.(this%m%b%g%nx==g%nx.and.abs(this%m%b%g%xmin-g%xmin)<eps.and.abs(this%m%b%g%xmax-g%xmax)<eps)) then
            stop "E: incompatible grids"
        end if
#elif(_DIM_==2)
        if (.not.(this%m%b%g%nx==g%nx.and.abs(this%m%b%g%xmin-g%xmin)<eps.and.abs(this%m%b%g%xmax-g%xmax)<eps &
             .and.this%m%b%g%ny==g%ny.and.abs(this%m%b%g%ymin-g%ymin)<eps.and.abs(this%m%b%g%ymax-g%ymax)<eps)) then
            stop "E: incompatible grids"
        end if
#elif(_DIM_==3)
        if (.not.(this%m%b%g%nx==g%nx.and.abs(this%m%b%g%xmin-g%xmin)<eps.and.abs(this%m%b%g%xmax-g%xmax)<eps &
             .and.this%m%b%g%ny==g%ny.and.abs(this%m%b%g%ymin-g%ymin)<eps.and.abs(this%m%b%g%ymax-g%ymax)<eps &
             .and.this%m%b%g%nz==g%nz.and.abs(this%m%b%g%zmin-g%zmin)<eps.and.abs(this%m%b%g%zmax-g%zmax)<eps)) then
            stop "E: incompatible grids"
        end if
#endif
        do k = 1, this%m%nc
            this%C(k)%p%u = 0.0_prec
#ifdef _REAL_        
            call hdf5_load_real_gridfun(this%m%b%g, this%C(k)%p%u, filename, trim(this%m%dset_names(k)))
#else
            call hdf5_load_complex_gridfun(this%m%b%g, this%C(k)%p%u, filename, &
                   trim(this%m%dset_names_real(k)),trim(this%m%dset_names_imag(k)))
#endif
        end do
#endif
    end subroutine S(load_wf_multicomponent_fourier)


    subroutine S(save_wf_multicomponent_fourier)(this, filename)
#ifdef _NO_HDF5_
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        character(len=*), intent(in) :: filename
        print *, "W: save not implemented"
#else
#ifdef _QUADPRECISION_
        use tssmq_hdf5
#else
        use tssm_hdf5
#endif        
        class(S(wf_multicomponent_fourier)), intent(inout) :: this
        character(len=*), intent(in) :: filename
        
        integer(HID_T) :: file_id  
        integer :: k

        call this%to_real_space()

        file_id = hdf5_open_gridfun(filename)
        do k = 1, this%m%nc
#ifdef _REAL_        
            call hdf5_write_real_gridfun(file_id, this%m%b%g, this%C(k)%p%u, trim(this%m%dset_names(k)))
#else
            call hdf5_write_complex_gridfun(file_id, this%m%b%g, this%C(k)%p%u, &
                   trim(this%m%dset_names_real(k)),trim(this%m%dset_names_imag(k)))
#endif
        end do
        call hdf5_close_gridfun(file_id)
#endif
    end subroutine S(save_wf_multicomponent_fourier)

#ifdef _QUADPRECISION_
end module S(tssmq_multicomponent_fourier)    
#else
end module S(tssm_multicomponent_fourier)    
#endif


