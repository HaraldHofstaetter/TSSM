
#if defined(_LAGUERRE_)

#ifdef _REAL_
#if (_DIM_==2)
#define SB(y) generalized_ ## y ## _real_2D
#define S(x) x ## _gen_lag_real_2D
#elif (_DIM_==3)
#define SB(y) generalized_ ## y ## _hermite_real_3D
#define S(x) x ## _gen_lag_herm_real_3D
#endif
#else
#if (_DIM_==2)
#define SB(y) generalized_ ## y ## _2D
#elif (_DIM_==3)
#define SB(y) generalized_ ## y ## _hermite_3D
#endif
#endif

#else

#if defined(_HERMITE_)
#ifdef _REAL_
#define S0(x,y)  x ## _hermite_real_ ## y ## D 
#else
#define S0(x,y)  x ## _hermite_ ## y ## D 
#endif
#else
#define _FOURIER_
#ifdef _REAL_
#define S0(x,y)  x ## _real_ ## y ## D 
#else
#define S0(x,y)  x ## _ ## y ## D 
#endif
#endif

#define S1(x,y) S0(x,y)
#define S(x) S1(x,_DIM_)

#ifdef _REAL_
#define SB0(x,y)  x ## _real_ ## y ## D 
#else
#define SB0(x,y)  x ## _ ## y ## D 
#endif
#define SB1(x,y) SB0(x,y)
#define SB(x) SB1(x,_DIM_)

#endif 

#ifdef _REAL_
#define _WAVE_FUNCTION_ real_wave_function
#define _COMPLEX_OR_REAL_ real
#else
#define _WAVE_FUNCTION_ wave_function
#define _COMPLEX_OR_REAL_ complex
#endif

#ifdef _QUADPRECISION_
#define fftw_destroy_plan fftwq_destroy_plan
#define fftw_free fftwq_free
#endif 



#ifdef _QUADPRECISION_
module S(tssmq_schroedinger) ! (Nonlinear) Schroedinger
    use tssmq_base
    use tssmq_fourier_common ! fftw_free
#if defined(_HERMITE_)
    use tssmq_hermite
#else
    use tssmq_fourier
#endif
#else
module S(tssm_schroedinger) ! (Nonlinear) Schroedinger
    use tssm_base
    use tssm_fourier_common ! fftw_free
#if defined(_HERMITE_)
    use tssm_hermite
#else
    use tssm_fourier
#endif
#endif
    implicit none

    private
    public :: S(schroedinger), S(wf_schroedinger)
#ifdef _QUADPRECISION_
    public :: wrapped_float128
#endif

#if defined(_HERMITE_)
    type, extends(SB(hermite)) :: S(schroedinger)
#else
    type, extends(SB(fourier)) :: S(schroedinger)
#endif
        real(kind=prec) :: hbar = 1.0_prec ! (default 1 for atomic units)
        real(kind=prec) :: mass = 1.0_prec ! (default 1 for atomic units)
#if defined(_HERMITE_)
        real(kind=prec) :: omega_x = 1.0_prec 
#if(_DIM_>=2)
        real(kind=prec) :: omega_y = 1.0_prec 
#endif
#if(_DIM_>=3)
        real(kind=prec) :: omega_z = 1.0_prec 
#endif
#endif
        real(kind=prec) :: cubic_coupling = 0.0_prec
#if(_DIM_==1)
        real(kind=prec), pointer :: V(:) => null()   ! for potential
        real(kind=prec), pointer :: V_t(:) => null() ! temporary storage for time-dependent potential 
#elif(_DIM_==2)
        real(kind=prec), pointer :: V(:,:) => null()
        real(kind=prec), pointer :: V_t(:,:) => null()
#elif(_DIM_==3)
        real(kind=prec), pointer :: V(:,:,:) => null()
        real(kind=prec), pointer :: V_t(:,:,:) => null()
#endif
        procedure(potential_t_interface), pointer, nopass :: potential_t => null()
        procedure(c_potential_t_interface), pointer, nopass :: c_potential_t => null()

        type(S(wf_schroedinger)), pointer :: tmp => null()

    contains
        procedure :: finalize => finalize_method
        procedure :: set_potential
#ifndef _REAL_        
        procedure :: set_potential_t
        procedure :: set_c_potential_t
#endif        
        procedure :: save_potential
        procedure :: load_potential
        procedure :: initialize_tmp
        procedure :: finalize_tmp
    end type S(schroedinger)

#ifdef _QUADPRECISION_
    type, bind(c) :: wrapped_float128
        real(kind=prec) :: v
    end type
#endif    

    abstract interface
#if(_DIM_==1)
        function potential_t_interface(x,t) 
#elif(_DIM_==2)
        function potential_t_interface(x,y,t) 
#elif(_DIM_==3)
        function potential_t_interface(x,y,z,t) 
#endif            
            import prec
#if(_DIM_==1)
            real(kind=prec), intent(in) :: x,t
#elif(_DIM_==2)
            real(kind=prec), intent(in) :: x,y,t
#elif(_DIM_==3)
            real(kind=prec), intent(in) :: x,y,z,t
#endif            
            real(kind=prec) :: potential_t_interface
        end function potential_t_interface

#if(_DIM_==1)
        function c_potential_t_interface(x,t) bind(c) 
#elif(_DIM_==2)
        function c_potential_t_interface(x,y,t) bind(c) 
#elif(_DIM_==3)
        function c_potential_t_interface(x,y,z,t) bind(c) 
#endif            
#ifdef _QUADPRECISION_
           import wrapped_float128
           type(wrapped_float128), value :: t
           type(wrapped_float128), value :: x
#if(_DIM_>=2)
           type(wrapped_float128), value :: y
#endif               
#if(_DIM_>=3)
           type(wrapped_float128), value :: z
#endif               
           type(wrapped_float128) :: c_potential_t_interface 
#else
            import prec
#if(_DIM_==1)
            real(kind=prec), value :: x,t
#elif(_DIM_==2)
            real(kind=prec), value :: x,y,t
#elif(_DIM_==3)
            real(kind=prec), value :: x,y,z,t
#endif            
            real(kind=prec) :: c_potential_t_interface
#endif            
        end function c_potential_t_interface
    end interface


    interface S(schroedinger)
        module procedure new_method 
    end interface S(schroedinger)


#if defined(_HERMITE_)
    type, extends(SB(wf_hermite)) :: S(wf_schroedinger)
#else
    type, extends(SB(wf_fourier)) :: S(wf_schroedinger)
#endif
    contains
        procedure :: clone
#ifndef _REAL_
        procedure :: propagate_B
#endif
        procedure :: imaginary_time_propagate_A
        procedure :: imaginary_time_propagate_B
        procedure :: add_apply_B
        procedure :: kinetic_energy
        procedure :: potential_energy
        procedure :: interaction_energy
        procedure :: observable
        procedure :: get_energy_expectation_deviation
        procedure :: get_realspace_observables
        procedure :: selfconsistent_nonlinear_step
        procedure :: extrapolation_imaginary_time_step
        procedure :: splitting_imaginary_time_step
        procedure :: compute_groundstate
        procedure :: print_local_imaginary_time_orders
    end type S(wf_schroedinger)

    interface S(wf_schroedinger)
        module procedure new_wf
    end interface S(wf_schroedinger)

contains

#if defined(_HERMITE_)

    function new_method( &
#if(_DIM_==1)
               nx, omega_x, &
#elif(_DIM_==2)
               nx, omega_x, ny, omega_y, &
#elif(_DIM_==3)
               nx, omega_x, ny, omega_y, nz, omega_z,  &
#endif
               hbar, mass, potential, &
               cubic_coupling) result(this)
        type(S(schroedinger)) :: this
        real(kind=prec), intent(in) :: omega_x 
        integer, intent(in) :: nx
#if(_DIM_>=2)
        real(kind=prec), intent(in) :: omega_y 
        integer, intent(in) :: ny
#endif
#if(_DIM_>=3)
        real(kind=prec), intent(in) :: omega_z 
        integer, intent(in) :: nz
#endif
        real(kind=prec), intent(in), optional :: hbar 
        real(kind=prec), intent(in), optional :: mass 
        real(kind=prec), external, optional   :: potential 
        real(kind=prec), intent(in), optional :: cubic_coupling 

        real(kind=prec) :: f 
        
        this%omega_x = omega_x
#if(_DIM_>=2)
        this%omega_y = omega_y
#endif
#if(_DIM_>=3)
        this%omega_z = omega_z
#endif

        f = sqrt(this%mass)/this%hbar
#if(_DIM_==1)
        this%SB(hermite) = SB(hermite)(nx, f*omega_x)
#elif(_DIM_==2)
        this%SB(hermite) = SB(hermite)(nx, f*omega_x, ny, f*omega_y)
#elif(_DIM_==3)
        this%SB(hermite) = SB(hermite)(nx, f*omega_x, ny, f*omega_y, nz, f*omega_z)
#endif
       if (present(hbar)) then
            this%hbar = hbar
        end if
        if (present(mass)) then
            this%mass = mass
        end if
        if (present(cubic_coupling)) then
            this%cubic_coupling = cubic_coupling
        end if

        if (present(potential)) then
            call this%set_potential(potential) 
        end if
    end function new_method


#else

    function new_method( &
#if(_DIM_==1)
               nx, xmin, xmax, &
#elif(_DIM_==2)
               nx, xmin, xmax, ny, ymin, ymax, &
#elif(_DIM_==3)
               nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax, &
#endif
               hbar, mass, potential, &
               cubic_coupling, boundary_conditions) result(this)
        type(S(schroedinger)) :: this
        real(kind=prec), intent(in) :: xmin 
        real(kind=prec), intent(in) :: xmax
        integer, intent(in) :: nx
#if(_DIM_>=2)
        real(kind=prec), intent(in) :: ymin 
        real(kind=prec), intent(in) :: ymax
        integer, intent(in) :: ny
#endif
#if(_DIM_>=3)
        real(kind=prec), intent(in) :: zmin 
        real(kind=prec), intent(in) :: zmax
        integer, intent(in) :: nz
#endif
        real(kind=prec), intent(in), optional :: hbar 
        real(kind=prec), intent(in), optional :: mass 
        real(kind=prec), external, optional   :: potential 
        real(kind=prec), intent(in), optional :: cubic_coupling 
        integer, intent(in), optional :: boundary_conditions

#if(_DIM_==1)
        this%SB(fourier) = SB(fourier)(nx, xmin, xmax, boundary_conditions)
#elif(_DIM_==2)
        this%SB(fourier) = SB(fourier)(nx, xmin, xmax, ny, ymin, ymax, boundary_conditions)
#elif(_DIM_==3)
        this%SB(fourier) = SB(fourier)(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax, boundary_conditions)
#endif
        if (present(hbar)) then
            this%hbar = hbar
        end if
        if (present(mass)) then
            this%mass = mass
        end if
        if (present(cubic_coupling)) then
            this%cubic_coupling = cubic_coupling
        end if
        if (present(potential)) then
            call this%set_potential(potential) 
        end if
    end function new_method

#endif


    subroutine finalize_method(this)
        class(S(schroedinger)), intent(inout) :: this

#if defined(_HERMITE_)
        call this%SB(hermite)%finalize()
#else
        call this%SB(fourier)%finalize()
#endif
        if(associated(this%V)) then 
            deallocate( this%V)
        end if 
        if(associated(this%V_t)) then 
            deallocate( this%V_t)
        end if 
        call this%finalize_tmp
    end subroutine finalize_method


    function new_wf(m) result(this)
        type(S(wf_schroedinger)) :: this
        class(S(schroedinger)), target, intent(inout) :: m
#if defined(_HERMITE_)
! TODO: nondefault values for hbar, omega
#ifdef _REAL_
        this%SB(wf_hermite) = SB(wf_hermite)(m, coefficient=m%hbar/m%mass)
#else
        this%SB(wf_hermite) = SB(wf_hermite)(m, coefficient=cmplx(0_prec, m%hbar/m%mass, kind=prec))
#endif
#else
#ifdef _REAL_
        this%SB(wf_fourier) = SB(wf_fourier)(m, coefficient=m%hbar/(2.0_prec*m%mass))
#else
        this%SB(wf_fourier) = SB(wf_fourier)(m, coefficient=cmplx(0.0_prec, m%hbar/(2.0_prec*m%mass), kind=prec))
#endif
#endif 
    end function new_wf


    function clone(this) 
        class(S(wf_schroedinger)), intent(inout) :: this
        class(_WAVE_FUNCTION_), pointer :: clone
        type(S(wf_schroedinger)), pointer :: p

        select type (m=>this%m); class is (S(schroedinger))
        allocate( p )
        p = S(wf_schroedinger)(m)
        clone => p
        end select
    end function clone



    subroutine set_potential(this, f)
        class(S(schroedinger)), intent(inout) :: this
        real(kind=prec), external :: f

        if (.not.associated(this%V)) then
        call this%g%allocate_real_gridfun(this%V)
        end if
        call this%g%set_real_gridfun(this%V, f)
    end subroutine set_potential

#ifndef _REAL_
    subroutine set_potential_t(this, f)
        class(S(schroedinger)), intent(inout) :: this
        real(kind=prec), external :: f
        if (.not.associated(this%V)) then
            call this%g%allocate_real_gridfun(this%V_t)
        end if
        this%potential_t => f
    end subroutine set_potential_t

    subroutine set_c_potential_t(this, f)
        class(S(schroedinger)), intent(inout) :: this
        interface 
#if(_DIM_==1)
           function f(x, t) bind(c)
#elif(_DIM_==2)
           function f(x, y, t) bind(c)
#elif(_DIM_==3)
           function f(x, y, z, t) bind(c)
#endif          
#ifdef _QUADPRECISION_
               import wrapped_float128
               type(wrapped_float128), value :: t
               type(wrapped_float128), value :: x
#if(_DIM_>=2)
               type(wrapped_float128), value :: y
#endif               
#if(_DIM_>=3)
               type(wrapped_float128), value :: z
#endif               
               type(wrapped_float128) :: f 
#else

               import prec
               real(kind=prec), value :: t
               real(kind=prec), value :: x
#if(_DIM_>=2)
               real(kind=prec), value :: y
#endif               
#if(_DIM_>=3)
               real(kind=prec), value :: z
#endif               
               real(kind=prec) :: f 
#endif               
           end function f
        end interface 
        this%c_potential_t => f
        call this%set_potential_t(eval_c_potential_t)
    contains
#if(_DIM_==1)    
        function eval_c_potential_t(x,t)
            real(kind=prec), intent(in) :: x,t
            real(kind=prec) :: eval_c_potential_t
#ifdef _QUADPRECISION_
            type(wrapped_float128) :: xx,tt,res
            xx%v = x
            tt%v = t
            res = this%c_potential_t(xx,tt)
            eval_c_potential_t = res%v
#else
            eval_c_potential_t = this%c_potential_t(x,t)
#endif            
        end function eval_c_potential_t
#elif(_DIM_==2)    
        function eval_c_potential_t(x,y,t)
            real(kind=prec), intent(in) :: x,y,t
            real(kind=prec) :: eval_c_potential_t
#ifdef _QUADPRECISION_
            type(wrapped_float128) :: xx,yy,tt,res
            xx%v = x
            yy%v = y
            tt%v = t
            res = this%c_potential_t(xx,yy,tt)
            eval_c_potential_t = res%v
#else
            eval_c_potential_t = this%c_potential_t(x,y,t)
#endif            
        end function eval_c_potential_t
#elif(_DIM_==3)    
        function eval_c_potential_t(x,y,z,t)
            real(kind=prec), intent(in) :: x,y,z,t
            real(kind=prec) :: eval_c_potential_t
#ifdef _QUADPRECISION_
            type(wrapped_float128) :: xx,yy,zz,tt,res
            xx%v = x
            yy%v = y
            zz%v = z
            tt%v = t
            res = this%c_potential_t(xx,yy,zz,tt)
            eval_c_potential_t = res%v
#else
            eval_c_potential_t = this%c_potential_t(x,y,z,t)
#endif        
        end function eval_c_potential_t
#endif        
    end subroutine set_c_potential_t
#endif


    subroutine save_potential(this, filename)
#ifdef _NO_HDF5_
        class(S(schroedinger)), intent(inout) :: this
        character(len=*), intent(in) :: filename
        print *, "W: save_potential not implemented"
#else
#ifdef _QUADPRECISION_
        use tssmq_hdf5
#else
        use tssm_hdf5
#endif        
        class(S(schroedinger)), intent(inout) :: this
        character(len=*), intent(in) :: filename
        call hdf5_save_real_gridfun(this%g, this%V, filename, "potential")
#if defined(_HERMITE_)
!TODO write grid attributes for hermite grid
#else
        call hdf5_write_grid_attributes(this%g, filename)
#endif
#endif
    end subroutine save_potential


!TODO slightly incompatible grids need not be allowed !
    subroutine load_potential(this, filename)
#ifdef _NO_HDF5_
        class(S(schroedinger)), intent(inout) :: this
        character(len=*), intent(in) :: filename
        print *, "W: load_potential not implemented"
#else
#ifdef _QUADPRECISION_
        use tssmq_hdf5
#else
        use tssm_hdf5
#endif        
        class(S(schroedinger)), intent(inout) :: this
        character(len=*), intent(in) :: filename
#ifdef _QUADPRECISION_
        real(kind=prec), parameter :: eps = epsilon(1.0_8)
#else        
        real(kind=prec), parameter :: eps = epsilon(1.0_prec)
#endif
        integer :: offset(3)
#if defined(_HERMITE_)
!TODO
#else
#if(_DIM_==1)
        type(grid_equidistant_1D) :: g 
 
#elif(_DIM_==2)
        type(grid_equidistant_2D) :: g 
#elif(_DIM_==3)
        type(grid_equidistant_3D) :: g 
#endif 
        call hdf5_read_grid_attributes(g, filename, "potential")
#if(_DIM_==1)
        if (.not.(abs(this%g%dx-g%dx)<eps)) then
            stop "E: incompatible grids (steps sizes do not match)"
        end if
        offset(1) = nint((g%xmin-this%g%xmin)/this%g%dx)
        this%g%xmin=g%xmin-offset(1)*this%g%dx
        this%g%xmax=this%g%xmin+this%g%nx*this%g%dx
#elif(_DIM_==2)
        if (.not.(abs(this%g%dx-g%dx)<eps.and.abs(this%g%dy-g%dy)<eps)) then
            stop "E: incompatible grids (steps sizes do not match)"
        end if
        offset(1) = nint((g%xmin-this%g%xmin)/this%g%dx)
        this%g%xmin=g%xmin-offset(1)*this%g%dx
        this%g%xmax=this%g%xmin+this%g%nx*this%g%dx
        offset(2) = nint((g%ymin-this%g%ymin)/this%g%dy)
        this%g%ymin=g%ymin-offset(2)*this%g%dy
        this%g%ymax=this%g%ymin+this%g%ny*this%g%dy
#elif(_DIM_==3)
        if (.not.(abs(this%g%dx-g%dx)<eps.and.abs(this%g%dy-g%dy)<eps.and.abs(this%g%dz-g%dz)<eps)) then
            stop "E: incompatible grids (steps sizes do not match)"
        end if
        offset(1) = nint((g%xmin-this%g%xmin)/this%g%dx)
        this%g%xmin=g%xmin-offset(1)*this%g%dx
        this%g%xmax=this%g%xmin+this%g%nx*this%g%dx
        offset(2) = nint((g%ymin-this%g%ymin)/this%g%dy)
        this%g%ymin=g%ymin-offset(2)*this%g%dy
        this%g%ymax=this%g%ymin+this%g%ny*this%g%dy
        offset(3) = nint((g%zmin-this%g%zmin)/this%g%dz)
        this%g%zmin=g%zmin-offset(3)*this%g%dz
        this%g%zmax=this%g%zmin+this%g%nz*this%g%dz
#endif
        if (.not.associated(this%V)) then
            call this%g%allocate_real_gridfun(this%V)
        end if
        this%V = 0.0_prec
        call hdf5_load_real_gridfun(this%g, this%V, filename, "potential", offset=offset)
#endif
#endif
    end subroutine load_potential




    subroutine imaginary_time_propagate_A(this, dt)
        class(S(wf_schroedinger)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt
#ifdef _REAL_
        call this%propagate_A(dt) 
#else
        call this%propagate_A(dt*cmplx(0.0_prec,-1.0_prec, kind=prec)) 
#endif
    end subroutine imaginary_time_propagate_A


    subroutine selfconsistent_nonlinear_step(this, dt, dt1, eps, max_iters) 
        class(S(wf_schroedinger)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt1
        real(kind=prec), intent(in), optional :: eps 
        integer, intent(in), optional :: max_iters

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
        integer :: iter, max_iters1
        real(kind=prec) :: N, N0, eps1
        

        eps1 =  100.0_prec*epsilon(1.0_prec)
        if (present(eps)) then 
           eps1 = eps
        end if   

        max_iters1 = 30
        if (present(max_iters)) then
            max_iters1 = max_iters
        end if    
        
        select type (m=>this%m); class is (S(schroedinger))

        call m%initialize_tmp
        call m%tmp%copy(this)

        N = m%tmp%norm()
!#ifndef _REAL_
!        if (aimag(dt1)>0.0_prec) then
!            this%u = exp((-aimag(dt1)/(N**2*m%hbar)*m%cubic_coupling*(0.0_prec, 1.0_prec)) &
!                         *(real(this%u, prec)**2+aimag(this%u)**2)) * this%u
!        end if                 
!#endif


        do iter = 1,max_iters1
            N0 = N
#ifndef _OPENMP
#ifdef _REAL_
            m%tmp%u = exp((-dt1*m%cubic_coupling/(N**2*m%hbar)) *  m%tmp%u**2 ) * this%u
#else            
            m%tmp%u = exp((-real(dt1, prec)*m%cubic_coupling/(N**2*m%hbar)) &
                           * ( real(m%tmp%u,prec)**2 + aimag(m%tmp%u)**2 ) ) * this%u
#endif         
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
            do j=1,n_threads
#if(_DIM_==1)
                u1 => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                u2 => m%tmp%u(lbound(m%tmp%u,1)+m%g%jj(j-1):lbound(m%tmp%u,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u1 => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                u2 => m%tmp%u(:,lbound(m%tmp%u,2)+m%g%jj(j-1):lbound(m%tmp%u,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u1 => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                u2 => m%tmp%u(:,:,lbound(m%tmp%u,3)+m%g%jj(j-1):lbound(m%tmp%u,3)+m%g%jj(j)-1)
#endif 
#ifdef _REAL_
                u2 = exp((-dt1*m%cubic_coupling/(N**2*m%hbar)) * u2**2 ) * u1
#else            
                u2 = exp((-real(dt1, prec)*m%cubic_coupling/(N**2*m%hbar)) &
                           * ( real(u2,prec)**2 + aimag(u2)**2 ) ) * u1
#endif         
            end do
!$OMP END PARALLEL DO
#endif    
            N = m%tmp%norm()
            if (abs(N-N0)<eps1) exit
        end do

#ifdef _REAL_
#ifndef _OPENMP
        this%u = exp((-dt*m%cubic_coupling/(N**2*m%hbar)) *  m%tmp%u**2 ) * this%u
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
            do j=1,n_threads
#if(_DIM_==1)
                u1 => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                u2 => m%tmp%u(lbound(m%tmp%u,1)+m%g%jj(j-1):lbound(m%tmp%u,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u1 => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                u2 => m%tmp%u(:,lbound(m%tmp%u,2)+m%g%jj(j-1):lbound(m%tmp%u,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u1 => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                u2 => m%tmp%u(:,:,lbound(m%tmp%u,3)+m%g%jj(j-1):lbound(m%tmp%u,3)+m%g%jj(j)-1)
#endif 
                u1 = exp((-dt*m%cubic_coupling/(N**2*m%hbar)) * u2**2 ) * u1
            end do
!$OMP END PARALLEL DO
#endif         
#else            
        if (aimag(dt)>0.0_prec) then
#ifndef _OPENMP
            this%u = exp(cmplx(0.0_prec, -aimag(dt1)/(N**2*m%hbar)*m%cubic_coupling,prec) &
                         *(real(m%tmp%u, prec)**2+aimag(m%tmp%u)**2)) * m%tmp%u
            this%u = exp((-real(dt-dt1, prec)*m%cubic_coupling/(N**2*m%hbar)) &
                      * ( real(this%u,prec)**2 + aimag(this%u)**2 ) ) * this%u
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
            do j=1,n_threads
#if(_DIM_==1)
                u1 => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                u2 => m%tmp%u(lbound(m%tmp%u,1)+m%g%jj(j-1):lbound(m%tmp%u,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u1 => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                u2 => m%tmp%u(:,lbound(m%tmp%u,2)+m%g%jj(j-1):lbound(m%tmp%u,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u1 => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                u2 => m%tmp%u(:,:,lbound(m%tmp%u,3)+m%g%jj(j-1):lbound(m%tmp%u,3)+m%g%jj(j)-1)
#endif 
                u1 = exp(cmplx(0.0_prec, -aimag(dt1)/(N**2*m%hbar)*m%cubic_coupling,prec) &
                         *(real(u2, prec)**2+aimag(u2)**2)) * u2
                u1 = exp((-real(dt-dt1, prec)*m%cubic_coupling/(N**2*m%hbar)) &
                      * ( real(u1,prec)**2 + aimag(u1)**2 ) ) * u1
            end do
!$OMP END PARALLEL DO
#endif
        else                 
#ifndef _OPENMP
            this%u = exp((-real(dt, prec)*m%cubic_coupling/(N**2*m%hbar)) &
                      * ( real(m%tmp%u,prec)**2 + aimag(m%tmp%u)**2 ) ) * this%u
#else
!$OMP PARALLEL DO PRIVATE(j, u1,u2) 
            do j=1,n_threads
#if(_DIM_==1)
                u1 => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                u2 => m%tmp%u(lbound(m%tmp%u,1)+m%g%jj(j-1):lbound(m%tmp%u,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u1 => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                u2 => m%tmp%u(:,lbound(m%tmp%u,2)+m%g%jj(j-1):lbound(m%tmp%u,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u1 => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                u2 => m%tmp%u(:,:,lbound(m%tmp%u,3)+m%g%jj(j-1):lbound(m%tmp%u,3)+m%g%jj(j)-1)
#endif 
                u1 = exp((-real(dt, prec)*m%cubic_coupling/(N**2*m%hbar)) &
                      * ( real(u2,prec)**2 + aimag(u2)**2 ) ) * u1
            end do
!$OMP END PARALLEL DO
#endif         
        end if                 
#endif            
!#ifndef _REAL_
!        N = this%norm()
!        if (aimag(dt-dt1)>0.0_prec) then
!            this%u = exp((-aimag(dt-dt1)/(N**2*m%hbar)*m%cubic_coupling*(0.0_prec, 1.0_prec)) &
!                         *(real(this%u, prec)**2+aimag(this%u)**2)) * this%u
!        end if                 
!#endif

        !if (abs(N-N0)>=eps1) then
        !     print *,  "W: iteration did not converge after", max_iters1, "steps"
        !end if     

        class default
           stop "E: wrong spectral method for schroedinger wave function"
        end select
    end subroutine selfconsistent_nonlinear_step


    subroutine imaginary_time_propagate_B(this, dt, method_for_B)
        class(S(wf_schroedinger)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt
        integer, intent(in), optional :: method_for_B
        integer :: method_for_B_1
#ifdef _OPENMP
#if(_DIM_==1)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:)
        real(kind=prec), pointer :: V(:)
#elif(_DIM_==2)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:)
        real(kind=prec), pointer :: V(:,:)
#elif(_DIM_==3)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:,:)
        real(kind=prec), pointer :: V(:,:,:)
#endif
        integer :: j
#endif        
        real(kind=prec) :: N

        if (dt==0.0_prec) then
            return
        endif
        
        select type (m=>this%m); class is (S(schroedinger))

        method_for_B_1 = 2
        if (present(method_for_B)) then
           method_for_B_1 = method_for_B
        end if   
        if (m%cubic_coupling==0.0_prec) then
           method_for_B_1 = 0
        end if   


        call this%to_real_space
            
        select case(method_for_B_1)
        case(0) ! use norm of left endpoint
           if (m%cubic_coupling/=0.0_prec) then
               N = this%norm()
#ifndef _OPENMP
#ifdef _REAL_               
               this%u = exp((-dt/m%hbar*m%cubic_coupling/N**2)*this%u**2) * this%u
#else               
               this%u = exp((-dt/m%hbar*m%cubic_coupling/N**2)*(real(this%u,prec)**2+aimag(this%u)**2)) * this%u
#endif               
#else
!$OMP PARALLEL DO PRIVATE(j, u) 
            do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
#endif 
#ifdef _REAL_               
               u = exp((-dt/m%hbar*m%cubic_coupling/N**2)*u**2) * u
#else               
               u = exp((-dt/m%hbar*m%cubic_coupling/N**2)*(real(u,prec)**2+aimag(u)**2)) * u
#endif               
            end do
!$OMP END PARALLEL DO
#endif 
           end if    

           if (associated(m%V)) then 
#ifndef _OPENMP
              this%u = exp((-dt/m%hbar)*m%V) * this%u
#else
!$OMP PARALLEL DO PRIVATE(j, u, V) 
            do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                V => m%V(lbound(m%V,1)+m%g%jj(j-1):lbound(m%V,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                V => m%V(:,lbound(m%V,2)+m%g%jj(j-1):lbound(m%V,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                V => m%V(:,:,lbound(m%V,3)+m%g%jj(j-1):lbound(m%V,3)+m%g%jj(j)-1)
#endif 
              u = exp((-dt/m%hbar)*V) * u
            end do
!$OMP END PARALLEL DO
#endif 
           end if

        case(1) ! use norm of right endpoint
           !Note: reverse order as in case(0)
           if (associated(m%V)) then
#ifndef _OPENMP
              this%u = exp((-dt/m%hbar)*m%V) * this%u
#else
!$OMP PARALLEL DO PRIVATE(j, u, V) 
            do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                V => m%V(lbound(m%V,1)+m%g%jj(j-1):lbound(m%V,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                V => m%V(:,lbound(m%V,2)+m%g%jj(j-1):lbound(m%V,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                V => m%V(:,:,lbound(m%V,3)+m%g%jj(j-1):lbound(m%V,3)+m%g%jj(j)-1)
#endif 
              u = exp((-dt/m%hbar)*V) * u
            end do
!$OMP END PARALLEL DO
#endif 
           end if

           if (m%cubic_coupling/=0.0_prec) then
               call this%selfconsistent_nonlinear_step(dt, dt)
           end if    

        case(2) ! use norm of midpoint

           if (associated(m%V)) then
#ifndef _OPENMP
              this%u = exp((-0.5_prec*dt/m%hbar)*m%V) * this%u
#else
!$OMP PARALLEL DO PRIVATE(j, u, V) 
            do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                V => m%V(lbound(m%V,1)+m%g%jj(j-1):lbound(m%V,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                V => m%V(:,lbound(m%V,2)+m%g%jj(j-1):lbound(m%V,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                V => m%V(:,:,lbound(m%V,3)+m%g%jj(j-1):lbound(m%V,3)+m%g%jj(j)-1)
#endif 
                u = exp((-0.5_prec*dt/m%hbar)*V) * u
            end do
!$OMP END PARALLEL DO
#endif 
           end if
           if (m%cubic_coupling/=0.0_prec) then
               call this%selfconsistent_nonlinear_step(dt, 0.5_prec*dt)
           end if    
           if (associated(m%V)) then
#ifndef _OPENMP
              this%u = exp((-0.5_prec*dt/m%hbar)*m%V) * this%u
#else
!$OMP PARALLEL DO PRIVATE(j, u, V) 
            do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                V => m%V(lbound(m%V,1)+m%g%jj(j-1):lbound(m%V,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                V => m%V(:,lbound(m%V,2)+m%g%jj(j-1):lbound(m%V,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                V => m%V(:,:,lbound(m%V,3)+m%g%jj(j-1):lbound(m%V,3)+m%g%jj(j)-1)
#endif 
                u = exp((-0.5_prec*dt/m%hbar)*V) * u
            end do
!$OMP END PARALLEL DO
#endif 
           end if

        case(3) ! use norm of left and right endpoint

           if (m%cubic_coupling/=0.0_prec) then
               N = this%norm()
#ifndef _OPENMP
#ifdef _REAL_               
               this%u = exp((-0.5_prec*dt/m%hbar*m%cubic_coupling/N**2)*this%u**2) * this%u
#else               
               this%u = exp((-0.5_prec*dt/m%hbar*m%cubic_coupling/N**2)*(real(this%u,prec)**2+aimag(this%u)**2)) * this%u
#endif   
#else
!$OMP PARALLEL DO PRIVATE(j, u) 
            do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
#endif 
#ifdef _REAL_               
               u = exp((-0.5_prec*dt/m%hbar*m%cubic_coupling/N**2)*u**2) * u
#else               
               u = exp((-0.5_prec*dt/m%hbar*m%cubic_coupling/N**2)*(real(u,prec)**2+aimag(u)**2)) * u
#endif   
            end do
!$OMP END PARALLEL DO
#endif 
           end if    

           if (associated(m%V)) then 
#ifndef _OPENMP
              this%u = exp((-dt/m%hbar)*m%V) * this%u
#else
!$OMP PARALLEL DO PRIVATE(j, u, V) 
            do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                V => m%V(lbound(m%V,1)+m%g%jj(j-1):lbound(m%V,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                V => m%V(:,lbound(m%V,2)+m%g%jj(j-1):lbound(m%V,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                V => m%V(:,:,lbound(m%V,3)+m%g%jj(j-1):lbound(m%V,3)+m%g%jj(j)-1)
#endif 
                u = exp((-dt/m%hbar)*V) * u
            end do
!$OMP END PARALLEL DO
#endif 
           end if

           if (m%cubic_coupling/=0.0_prec) then
               call this%selfconsistent_nonlinear_step(0.5_prec*dt, 0.5_prec*dt)
           end if    

        case default
            stop "E: method for B not implemented"
        end select


        class default
           stop "E: wrong spectral method for schroedinger wave function"
        end select
    end subroutine imaginary_time_propagate_B



#ifndef _REAL_
    
    subroutine propagate_B(this, dt)
        class(S(wf_schroedinger)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt
        complex(kind=prec) :: f
#ifdef _OPENMP
#if(_DIM_==1)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:)
        real(kind=prec), pointer :: V(:)
#elif(_DIM_==2)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:)
        real(kind=prec), pointer :: V(:,:)
#elif(_DIM_==3)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:,:)
        real(kind=prec), pointer :: V(:,:,:)
#endif
        integer :: j
#endif        
        if (dt==0.0_prec) then
            return
        endif
        
        select type (m=>this%m); class is (S(schroedinger))
        call this%to_real_space

        f = -dt/m%hbar*(0.0_prec, 1.0_prec)

        if (associated(m%V)) then
!!!! CHECK Speicherzugriffsfehler
#ifndef _OPENMP
           this%u = exp(f*m%V) * this%u
#else
!$OMP PARALLEL DO PRIVATE(j, u, V) 
            do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                V => m%V(lbound(m%V,1)+m%g%jj(j-1):lbound(m%V,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                V => m%V(:,lbound(m%V,2)+m%g%jj(j-1):lbound(m%V,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                V => m%V(:,:,lbound(m%V,3)+m%g%jj(j-1):lbound(m%V,3)+m%g%jj(j)-1)
#endif 
                u = exp(f*V) * u
            end do
!$OMP END PARALLEL DO
#endif 
        end if

        if (associated(m%potential_t)) then
            call m%g%set_t_real_gridfun(m%V_t, m%potential_t, this%time)
#ifndef _OPENMP
            this%u = exp(f*m%V_t) * this%u
#else
!$OMP PARALLEL DO PRIVATE(j, u, V) 
            do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                V => m%V_t(lbound(m%V_t,1)+m%g%jj(j-1):lbound(m%V_t,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                V => m%V_t(:,lbound(m%V_t,2)+m%g%jj(j-1):lbound(m%V_t,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                V => m%V_t(:,:,lbound(m%V_t,3)+m%g%jj(j-1):lbound(m%V_t,3)+m%g%jj(j)-1)
#endif 
                u = exp(f*V) * u
            end do
!$OMP END PARALLEL DO
#endif 
        end if

        if (m%cubic_coupling/=0.0_prec) then
#ifndef _OPENMP
           this%u = exp(m%cubic_coupling*f*(real(this%u, prec)**2+aimag(this%u)**2)) * this%u
#else
!$OMP PARALLEL DO PRIVATE(j, u) 
            do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
#endif 
                u = exp(m%cubic_coupling*f*(real(u, prec)**2+aimag(u)**2)) * u
            end do
!$OMP END PARALLEL DO
#endif 
        end if    

        class default
           stop "E: wrong spectral method for schroedinger wave function"
        end select
    end subroutine propagate_B

#endif


    subroutine initialize_tmp(this)
        class(S(schroedinger)), intent(inout) :: this
        if (associated(this%tmp)) return 
        allocate( this%tmp )
        this%tmp = S(wf_schroedinger)(this)
    end subroutine initialize_tmp


    subroutine finalize_tmp(this)
        use, intrinsic :: iso_c_binding, only: c_loc
        class(S(schroedinger)), intent(inout) :: this
        if (.not.associated(this%tmp)) return
#if defined(_FOURIER_)
        call fftw_destroy_plan(this%tmp%plan_forward)
        call fftw_destroy_plan(this%tmp%plan_backward)
#endif
        !call fftw_free(c_loc(this%tmp%up))
        call fftw_free(c_loc(this%tmp%up(1)))
        deallocate( this%tmp )
        this%tmp => null()
    end subroutine finalize_tmp


    subroutine add_apply_B(this, wf, coefficient)
        class(S(wf_schroedinger)), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: wf 
        _COMPLEX_OR_REAL_(kind=prec), intent(in), optional :: coefficient
        _COMPLEX_OR_REAL_(kind=prec) :: C 
#ifdef _OPENMP
#if(_DIM_==1)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u1(:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u2(:)
        real(kind=prec), pointer :: V(:)
#elif(_DIM_==2)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u1(:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u2(:,:)
        real(kind=prec), pointer :: V(:,:)
#elif(_DIM_==3)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u1(:,:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u2(:,:,:)
        real(kind=prec), pointer :: V(:,:,:)
#endif
        integer :: j
#endif        


        select type (wf)
        class is (S(wf_schroedinger))
        
        if (.not.associated(wf%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    

        C = 1.0_prec
        if (present(coefficient)) then
            C = coefficient
        end if    

        select type (m=>this%m); class is (S(schroedinger))
        call this%to_real_space
        call wf%to_real_space

        if (associated(m%V)) then
#ifndef _OPENMP
#ifdef _REAL_
           wf%u = wf%u + (-C/m%hbar) * m%V * this%u
#else
           wf%u = wf%u + (-C/m%hbar*(0.0_prec, 1.0_prec))* m%V * this%u
#endif 
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2, V) 
           do j=1,n_threads
#if(_DIM_==1)
                u1 => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                u2 => wf%u(lbound(wf%u,1)+m%g%jj(j-1):lbound(wf%u,1)+m%g%jj(j)-1)
                V => m%V(lbound(m%V,1)+m%g%jj(j-1):lbound(m%V,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u1 => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                u2 => wf%u(:,lbound(wf%u,2)+m%g%jj(j-1):lbound(wf%u,2)+m%g%jj(j)-1)
                V => m%V(:,lbound(m%V,2)+m%g%jj(j-1):lbound(m%V,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u1 => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                u2 => wf%u(:,:,lbound(wf%u,3)+m%g%jj(j-1):lbound(wf%u,3)+m%g%jj(j)-1)
                V => m%V(:,:,lbound(m%V,3)+m%g%jj(j-1):lbound(m%V,3)+m%g%jj(j)-1)
#endif 
#ifdef _REAL_
                u2 = u2 + (-C/m%hbar) * V * u1
#else
                u2 = u2 + (-C/m%hbar*(0.0_prec, 1.0_prec))* V * u1
#endif 
           end do
!$OMP END PARALLEL DO
#endif 
        end if

        if (m%cubic_coupling/=0_prec) then
#ifndef _OPENMP
#ifdef _REAL_
           wf%u = wf%u + ((-C/m%hbar*m%cubic_coupling) * this%u**3) 
#else
           wf%u = wf%u + ((-C/m%hbar*m%cubic_coupling*(0.0_prec, 1.0_prec)) &
                         *(real(this%u, prec)**2+aimag(this%u)**2)) * this%u
#endif 
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
           do j=1,n_threads
#if(_DIM_==1)
                u1 => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                u2 => wf%u(lbound(wf%u,1)+m%g%jj(j-1):lbound(wf%u,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u1 => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                u2 => wf%u(:,lbound(wf%u,2)+m%g%jj(j-1):lbound(wf%u,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u1 => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                u2 => wf%u(:,:,lbound(wf%u,3)+m%g%jj(j-1):lbound(wf%u,3)+m%g%jj(j)-1)
#endif 
#ifdef _REAL_
                u2 = u2 + ((-C/m%hbar*m%cubic_coupling) * u1**3) 
#else
                u2 = u2 + ((-C/m%hbar*m%cubic_coupling*(0.0_prec, 1.0_prec)) &
                          *(real(u1, prec)**2+aimag(u1)**2)) * u1
#endif 
           end do
!$OMP END PARALLEL DO
#endif 
        end if   

        end select

        class default
           stop "E: wrong spectral method for schroedinger wave function"
        end select
    end subroutine add_apply_B
 


    function kinetic_energy(this) result(E_kin)
        class(S(wf_schroedinger)), intent(inout) :: this
        real(kind=prec) :: E_kin
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
        real(kind=prec), pointer :: ev(:)
        integer :: j
#endif        

        select type (m=>this%m); class is (S(schroedinger))

        call this%to_frequency_space
#ifdef _REAL_
#if defined(_HERMITE_)
#if(_DIM_==1)
#ifndef _OPENMP
        h = sum(this%m%eigenvalues1 &
            * this%uf**2 )
#else
        h = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, uf, ev) REDUCTION(+:h)
        do j=1,n_threads
            uf => this%uf(lbound(this%uf,1)+this%m%jf(j-1):lbound(this%uf,1)+this%m%jf(j)-1)
            ev => this%m%eigenvalues1(lbound(this%m%eigenvalues1,1)+this%m%jf(j-1):&
                                         lbound(this%m%eigenvalues1,1)+this%m%jf(j)-1)
            h = h +  sum(ev &
                     * uf**2 )
        end do
!$OMP END PARALLEL DO 
#endif
#elif(_DIM_==2)
#ifndef _OPENMP
        h = sum( (spread(this%m%eigenvalues1,2, this%m%nf2max-this%m%nf2min+1) &
                + spread(this%m%eigenvalues2,1, this%m%nf1max-this%m%nf1min+1)) &
            * this%uf**2 )
#else
        h = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, uf, ev) REDUCTION(+:h)
        do j=1,n_threads
            uf => this%uf(:,lbound(this%uf,2)+this%m%jf(j-1):lbound(this%uf,2)+this%m%jf(j)-1)
            ev => this%m%eigenvalues2(lbound(this%m%eigenvalues2,1)+this%m%jf(j-1):&
                                         lbound(this%m%eigenvalues2,1)+this%m%jf(j)-1)
            h = h + sum( (spread(this%m%eigenvalues1,2, this%m%jf(j)-this%m%jf(j-1)) &
                       + spread(ev,1, this%m%nf1max-this%m%nf1min+1)) &
                       * uf**2 )
        end do
!$OMP END PARALLEL DO 
#endif
#elif(_DIM_==3)
#ifndef _OPENMP
        h = sum((spread(spread(this%m%eigenvalues1, 2,this%m%nf2max-this%m%nf2min+1), 3,this%m%nf3max-this%m%nf3min+1) &
               + spread(spread(this%m%eigenvalues2, 1,this%m%nf1max-this%m%nf1min+1), 3,this%m%nf3max-this%m%nf3min+1) &
               + spread(spread(this%m%eigenvalues3, 1,this%m%nf1max-this%m%nf1min+1), 2,this%m%nf2max-this%m%nf2min+1)) &
            * this%uf**2 )
#else
        h = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, uf, ev) REDUCTION(+:h)
        do j=1,n_threads
            uf => this%uf(:,:,lbound(this%uf,3)+this%m%jf(j-1):lbound(this%uf,3)+this%m%jf(j)-1)
            ev => this%m%eigenvalues3(lbound(this%m%eigenvalues3,1)+this%m%jf(j-1):&
                                         lbound(this%m%eigenvalues3,1)+this%m%jf(j)-1)
            h = h + sum((spread(spread(this%m%eigenvalues1, 2,this%m%nf2max-this%m%nf2min+1), &
                                       3,this%m%jf(j)-this%m%jf(j-1)) &
                  + spread(spread(this%m%eigenvalues2, 1,this%m%nf1max-this%m%nf1min+1), &
                                       3,this%m%jf(j)-this%m%jf(j-1)) &
                  + spread(spread(ev, 1,this%m%nf1max-this%m%nf1min+1), 2,this%m%nf2max-this%m%nf2min+1)) &
                  * uf**2 )
        end do
!$OMP END PARALLEL DO 
#endif
#endif

#else
        select case(this%m%boundary_conditions)
        case(periodic) 
        this%m%eigenvalues1(this%m%nf1max) = 0.5_prec * this%m%eigenvalues1(this%m%nf1max)
#if(_DIM_==1)
#ifndef _OPENMP
        h = 2.0_prec*sum(this%m%eigenvalues1 * (real(this%ufc,prec)**2 + aimag(this%ufc)**2) ) 
#else
        h = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, ufc, ev) REDUCTION(+:h)
        do j=1,n_threads
            ufc => this%ufc(lbound(this%ufc,1)+this%m%jf(j-1):lbound(this%ufc,1)+this%m%jf(j)-1)
            ev => this%m%eigenvalues1(lbound(this%m%eigenvalues1,1)+this%m%jf(j-1):&
                                         lbound(this%m%eigenvalues1,1)+this%m%jf(j)-1)
            h = h + 2.0_prec*sum(ev * (real(ufc,prec)**2 + aimag(ufc)**2) ) 
        end do
!$OMP END PARALLEL DO 
#endif 
#elif(_DIM_==2)
#ifndef _OPENMP
        h = 2.0_prec*sum( (spread(this%m%eigenvalues1,2, this%m%nf2max-this%m%nf2min+1) &
                + spread(this%m%eigenvalues2,1, this%m%nf1max-this%m%nf1min+1)) &
            * (real(this%ufc,prec)**2 + aimag(this%ufc)**2) ) 
#else
        h = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, ufc, ev) REDUCTION(+:h)
        do j=1,n_threads
            ufc => this%ufc(:,lbound(this%ufc,2)+this%m%jf(j-1):lbound(this%ufc,2)+this%m%jf(j)-1)
            ev => this%m%eigenvalues2(lbound(this%m%eigenvalues2,1)+this%m%jf(j-1):&
                                         lbound(this%m%eigenvalues2,1)+this%m%jf(j)-1)
        h = h + 2.0_prec*sum( (spread(this%m%eigenvalues1,2, this%m%jf(j)-this%m%jf(j-1) ) &
                             + spread(ev,1, this%m%nf1max-this%m%nf1min+1)) &
                             * (real(ufc,prec)**2 + aimag(ufc)**2) ) 
        end do
!$OMP END PARALLEL DO 
#endif 
#elif(_DIM_==3)
#ifndef _OPENMP
        h = 2.0_prec &
          *sum((spread(spread(this%m%eigenvalues1, 2,this%m%nf2max-this%m%nf2min+1), 3,this%m%nf3max-this%m%nf3min+1) &
             + spread(spread(this%m%eigenvalues2, 1,this%m%nf1max-this%m%nf1min+1), 3,this%m%nf3max-this%m%nf3min+1) &
             + spread(spread(this%m%eigenvalues3, 1,this%m%nf1max-this%m%nf1min+1), 2,this%m%nf2max-this%m%nf2min+1))& 
            * (real(this%ufc,prec)**2 + aimag(this%ufc)**2) ) 
#else
        h = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, ufc, ev) REDUCTION(+:h)
        do j=1,n_threads
            ufc => this%ufc(:,:,lbound(this%ufc,3)+this%m%jf(j-1):lbound(this%ufc,3)+this%m%jf(j)-1)
            ev => this%m%eigenvalues3(lbound(this%m%eigenvalues3,1)+this%m%jf(j-1):&
                                         lbound(this%m%eigenvalues3,1)+this%m%jf(j)-1)
            h = h +  2.0_prec &
          *sum((spread(spread(this%m%eigenvalues1, 2,this%m%nf2max-this%m%nf2min+1), 3,  this%m%jf(j)-this%m%jf(j-1)) &
             + spread(spread(this%m%eigenvalues2, 1,this%m%nf1max-this%m%nf1min+1), 3, this%m%jf(j)-this%m%jf(j-1) ) &
             + spread(spread(ev, 1,this%m%nf1max-this%m%nf1min+1), 2,this%m%nf2max-this%m%nf2min+1))& 
            * (real(ufc,prec)**2 + aimag(ufc)**2) ) 
        end do
!$OMP END PARALLEL DO 
#endif 
#endif
        this%m%eigenvalues1(this%m%nf1max) = 2.0_prec * this%m%eigenvalues1(this%m%nf1max)
        case(dirichlet, neumann)
#if(_DIM_==1)
#ifndef _OPENMP
        h = sum(this%m%eigenvalues1 &
            * this%uf**2 )
#else
        h = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, uf, ev) REDUCTION(+:h)
        do j=1,n_threads
            uf => this%uf(lbound(this%uf,1)+this%m%jf(j-1):lbound(this%uf,1)+this%m%jf(j)-1)
            ev => this%m%eigenvalues1(lbound(this%m%eigenvalues1,1)+this%m%jf(j-1):&
                                         lbound(this%m%eigenvalues1,1)+this%m%jf(j)-1)
            h = h +  sum(ev &
                     * uf**2 )
        end do
!$OMP END PARALLEL DO 
#endif
#elif(_DIM_==2)
#ifndef _OPENMP
        h = sum( (spread(this%m%eigenvalues1,2, this%m%nf2max-this%m%nf2min+1) &
                + spread(this%m%eigenvalues2,1, this%m%nf1max-this%m%nf1min+1)) &
            * this%uf**2 )
#else
        h = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, uf, ev) REDUCTION(+:h)
        do j=1,n_threads
            uf => this%uf(:,lbound(this%uf,2)+this%m%jf(j-1):lbound(this%uf,2)+this%m%jf(j)-1)
            ev => this%m%eigenvalues2(lbound(this%m%eigenvalues2,1)+this%m%jf(j-1):&
                                         lbound(this%m%eigenvalues2,1)+this%m%jf(j)-1)
            h = h + sum( (spread(this%m%eigenvalues1,2, this%m%jf(j)-this%m%jf(j-1)) &
                       + spread(ev,1, this%m%nf1max-this%m%nf1min+1)) &
                       * uf**2 )
        end do
!$OMP END PARALLEL DO 
#endif
#elif(_DIM_==3)
#ifndef _OPENMP
        h = sum((spread(spread(this%m%eigenvalues1, 2,this%m%nf2max-this%m%nf2min+1), 3,this%m%nf3max-this%m%nf3min+1) &
               + spread(spread(this%m%eigenvalues2, 1,this%m%nf1max-this%m%nf1min+1), 3,this%m%nf3max-this%m%nf3min+1) &
               + spread(spread(this%m%eigenvalues3, 1,this%m%nf1max-this%m%nf1min+1), 2,this%m%nf2max-this%m%nf2min+1)) &
            * this%uf**2 )
#else
        h = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, uf, ev) REDUCTION(+:h)
        do j=1,n_threads
            uf => this%uf(:,:,lbound(this%uf,3)+this%m%jf(j-1):lbound(this%uf,3)+this%m%jf(j)-1)
            ev => this%m%eigenvalues3(lbound(this%m%eigenvalues3,1)+this%m%jf(j-1):&
                                         lbound(this%m%eigenvalues3,1)+this%m%jf(j)-1)
            h = h + sum((spread(spread(this%m%eigenvalues1, 2,this%m%nf2max-this%m%nf2min+1), &
                                       3,this%m%jf(j)-this%m%jf(j-1)) &
                  + spread(spread(this%m%eigenvalues2, 1,this%m%nf1max-this%m%nf1min+1), &
                                       3,this%m%jf(j)-this%m%jf(j-1)) &
                  + spread(spread(ev, 1,this%m%nf1max-this%m%nf1min+1), 2,this%m%nf2max-this%m%nf2min+1)) &
                  * uf**2 )
        end do
!$OMP END PARALLEL DO 
#endif
#endif
        end select
#endif
#else
#if(_DIM_==1)
#ifndef _OPENMP
        h = sum(this%m%eigenvalues1 &
            * (real(this%uf,prec)**2 + aimag(this%uf)**2) )
#else
        h = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, uf, ev) REDUCTION(+:h)
        do j=1,n_threads
            uf => this%uf(lbound(this%uf,1)+this%m%jf(j-1):lbound(this%uf,1)+this%m%jf(j)-1)
            ev => this%m%eigenvalues1(lbound(this%m%eigenvalues1,1)+this%m%jf(j-1):&
                                         lbound(this%m%eigenvalues1,1)+this%m%jf(j)-1)
            h = h +  sum(ev &
                  * (real(uf,prec)**2 + aimag(uf)**2) )
        end do
!$OMP END PARALLEL DO 
#endif
#elif(_DIM_==2)
#ifndef _OPENMP
        h = sum( (spread(this%m%eigenvalues1,2, this%m%nf2max-this%m%nf2min+1) &
                + spread(this%m%eigenvalues2,1, this%m%nf1max-this%m%nf1min+1)) &
                * (real(this%uf,prec)**2 + aimag(this%uf)**2) )
#else
        h = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, uf, ev) REDUCTION(+:h)
        do j=1,n_threads
            uf => this%uf(:,lbound(this%uf,2)+this%m%jf(j-1):lbound(this%uf,2)+this%m%jf(j)-1)
            ev => this%m%eigenvalues2(lbound(this%m%eigenvalues2,1)+this%m%jf(j-1):&
                                         lbound(this%m%eigenvalues2,1)+this%m%jf(j)-1)
            h = h + sum( (spread(this%m%eigenvalues1,2, this%m%jf(j)-this%m%jf(j-1)) &
                       + spread(ev,1, this%m%nf1max-this%m%nf1min+1)) &
                       * (real(uf,prec)**2 + aimag(uf)**2) )
        end do
!$OMP END PARALLEL DO 
#endif
#elif(_DIM_==3)
#ifndef _OPENMP
        h = sum((spread(spread(this%m%eigenvalues1, 2,this%m%nf2max-this%m%nf2min+1), 3,this%m%nf3max-this%m%nf3min+1) &
               + spread(spread(this%m%eigenvalues2, 1,this%m%nf1max-this%m%nf1min+1), 3,this%m%nf3max-this%m%nf3min+1) &
               + spread(spread(this%m%eigenvalues3, 1,this%m%nf1max-this%m%nf1min+1), 2,this%m%nf2max-this%m%nf2min+1)) &
               * (real(this%uf,prec)**2 + aimag(this%uf)**2) )
#else
        h = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, uf, ev) REDUCTION(+:h)
        do j=1,n_threads
            uf => this%uf(:,:,lbound(this%uf,3)+this%m%jf(j-1):lbound(this%uf,3)+this%m%jf(j)-1)
            ev => this%m%eigenvalues3(lbound(this%m%eigenvalues3,1)+this%m%jf(j-1):&
                                         lbound(this%m%eigenvalues3,1)+this%m%jf(j)-1)
            h = h + sum((spread(spread(this%m%eigenvalues1, 2,this%m%nf2max-this%m%nf2min+1), &
                                       3,this%m%jf(j)-this%m%jf(j-1)) &
                  + spread(spread(this%m%eigenvalues2, 1,this%m%nf1max-this%m%nf1min+1), &
                                       3,this%m%jf(j)-this%m%jf(j-1)) &
                  + spread(spread(ev, 1,this%m%nf1max-this%m%nf1min+1), 2,this%m%nf2max-this%m%nf2min+1)) &
                  * (real(uf,prec)**2 + aimag(uf)**2) )
        end do
!$OMP END PARALLEL DO 
#endif
#endif
#endif
#ifdef _MPI_
        call MPI_Reduce(h, h1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            h = h1
#endif
#if defined(_HERMITE_)
            E_kin =  -m%hbar**2/m%mass * h 
#else
#if(_DIM_==1)
            E_kin =  -0.5_prec*m%hbar**2/m%mass * h * m%g%dx/(m%g%nx)                                        
#elif(_DIM_==2)
            E_kin =  -0.5_prec*m%hbar**2/m%mass * h * m%g%dx*m%g%dy/(m%g%nx*m%g%ny)                                         
#elif(_DIM_==3)
            E_kin =  -0.5_prec*m%hbar**2/m%mass * h * m%g%dx*m%g%dy*m%g%dz/(m%g%nx*m%g%ny*m%g%nz) 
#endif
            select case(this%m%boundary_conditions)
            case(periodic)
            case(dirichlet, neumann)
                E_kin = E_kin/(2**_DIM_)
            end select
#endif
#ifdef _MPI_
        end if 
        call MPI_Bcast(E_kin, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
        class default
           stop "E: wrong spectral method for schroedinger wave function"
        end select
    end function kinetic_energy



    function potential_energy(this) result(E_pot)
        class(S(wf_schroedinger)), intent(inout) :: this
        real(kind=prec) :: E_pot
        real(kind=prec) :: dV 
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: h
#endif       
        integer :: n1, n2, n3
        integer :: nv1, nv2, nv3
#ifdef _OPENMP
#if(_DIM_==1)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:)
        real(kind=prec), pointer :: V(:)
#elif(_DIM_==2)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:)
        real(kind=prec), pointer :: V(:,:)
#elif(_DIM_==3)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:,:)
        real(kind=prec), pointer :: V(:,:,:)
#endif
#if defined(_HERMITE_)
        real(kind=prec), pointer :: w(:)
#endif        
        integer :: j
#endif        


        select type (m=>this%m); class is (S(schroedinger))

        if (associated(m%V)) then
           call this%to_real_space

#if defined(_HERMITE_)
#ifndef _OPENMP
#if(_DIM_==1)
          E_pot = sum(m%g%weights_x &
#elif(_DIM_==2)
          E_pot = sum(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1) &
                      *spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1) &
#elif(_DIM_==3)
          E_pot = sum(spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), 3, m%g%n3max-m%g%n3min+1) &
                     *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), 3, m%g%n3max-m%g%n3min+1) &
                     *spread(spread(m%g%weights_z, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) &
#endif
#ifdef _REAL_
                     *m%V *(this%u**2) )
#else
                     *m%V *(real(this%u,kind=prec)**2 +aimag(this%u)**2) )
#endif 
#else
          E_pot = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, V, w) REDUCTION(+:E_pot)
          do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                V => m%V(lbound(m%V,1)+m%g%jj(j-1):lbound(m%V,1)+m%g%jj(j)-1)
                w => m%g%weights_x(lbound(m%g%weights_x,1)+m%g%jj(j-1):lbound(m%g%weights_x,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                V => m%V(:,lbound(m%V,2)+m%g%jj(j-1):lbound(m%V,2)+m%g%jj(j)-1)
                w => m%g%weights_y(lbound(m%g%weights_y,1)+m%g%jj(j-1):lbound(m%g%weights_y,1)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                V => m%V(:,:,lbound(m%V,3)+m%g%jj(j-1):lbound(m%V,3)+m%g%jj(j)-1)
                w => m%g%weights_z(lbound(m%g%weights_z,1)+m%g%jj(j-1):lbound(m%g%weights_z,1)+m%g%jj(j)-1)
#endif 
#if(_DIM_==1)
                E_pot = E_pot + sum(w &
#elif(_DIM_==2)
                E_pot = E_pot + sum(spread(m%g%weights_x, 2, m%g%jj(j)-m%g%jj(j-1)) &
                      *spread(w, 1, m%g%n1max-m%g%n1min+1) &
#elif(_DIM_==3)
                E_pot = E_pot + sum(spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), 3, m%g%jj(j)-m%g%jj(j-1)) &
                     *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), 3, m%g%jj(j)-m%g%jj(j-1)) &
                     *spread(spread(w, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) &
#endif
#ifdef _REAL_
                     *V *(u**2) )
#else
                     *V *(real(u,kind=prec)**2 +aimag(u)**2) )
#endif 
          end do
!$OMP END PARALLEL DO
#endif 
#else
#if(_DIM_==1)
            dV = m%g%dx 
#elif(_DIM_==2)
            dV = m%g%dx * m%g%dy 
#elif(_DIM_==3)
            dV = m%g%dx * m%g%dy * m%g%dz
#endif
            select case(m%boundary_conditions)
            case(periodic, dirichlet)
#ifndef _OPENMP
#ifdef _REAL_
                E_pot = sum(m%V * this%u**2) * dV 
#else
                E_pot = sum(m%V * (real(this%u, prec)**2+aimag(this%u)**2)) * dV 
#endif 
#else
          E_pot = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, V) REDUCTION(+:E_pot)
          do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                V => m%V(lbound(m%V,1)+m%g%jj(j-1):lbound(m%V,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                V => m%V(:,lbound(m%V,2)+m%g%jj(j-1):lbound(m%V,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                V => m%V(:,:,lbound(m%V,3)+m%g%jj(j-1):lbound(m%V,3)+m%g%jj(j)-1)
#endif
#ifdef _REAL_
                E_pot = E_pot + sum(V * u**2) * dV 
#else
                E_pot = E_pot + sum(V * (real(u, prec)**2+aimag(u)**2)) * dV 
#endif 
          end do
!$OMP END PARALLEL DO
#endif 
            case(neumann) ! for summation exlude first index in each dimension 
#if(_DIM_==1)
                n1 = 1;  if (m%g%n1min==0) n1 = 2
                nv1 = lbound(m%V,1) + n1 - 1
#ifndef _OPENMP
#ifdef _REAL_
                E_pot = sum(m%V(nv1:) * this%u(n1:)**2) * dV 
#else
                E_pot = sum(m%V(nv1:) * (real(this%u(n1:), prec)**2+aimag(this%u(n1:))**2)) * dV 
#endif 
#else
                E_pot = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, V) REDUCTION(+:E_pot)
                do j=1,n_threads
                      u => this%u(max(n1,lbound(this%u,1)+m%g%jj(j-1)):lbound(this%u,1)+m%g%jj(j)-1)
                      V => m%V(max(nv1,lbound(m%V,1)+m%g%jj(j-1)):lbound(m%V,1)+m%g%jj(j)-1)
#ifdef _REAL_
                      E_pot = E_pot + sum(V * u**2) * dV 
#else
                      E_pot = E_pot + sum(V * (real(u, prec)**2+aimag(u)**2)) * dV 
#endif 
                end do
!$OMP END PARALLEL DO
#endif 
#elif(_DIM_==2)
                n1 = 1;  if (m%g%n1min==0) n1 = 2
                n2 = 1;  if (m%g%n2min==0) n2 = 2
                nv1 = lbound(m%V,1) + n1 - 1
                nv2 = lbound(m%V,2) + n2 - 1
#ifndef _OPENMP
#ifdef _REAL_
                E_pot = sum(m%V(nv1:,nv2:) * this%u(n1:,n2:)**2) * dV 
#else
                E_pot = sum(m%V(nv1:,nv2:) * (real(this%u(n1:,n2:), prec)**2 &
                                            +aimag(this%u(n1:,n2:))**2)) * dV 
#endif 
#else
                E_pot = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, V) REDUCTION(+:E_pot)
                do j=1,n_threads
                      u => this%u(n1:,max(n2,lbound(this%u,2)+m%g%jj(j-1)):lbound(this%u,2)+m%g%jj(j)-1)
                      V => m%V(nv1:,max(nv2,lbound(m%V,2)+m%g%jj(j-1)):lbound(m%V,2)+m%g%jj(j)-1)
#ifdef _REAL_
                      E_pot = E_pot + sum(V * u**2) * dV 
#else
                      E_pot = E_pot + sum(V * (real(u, prec)**2+aimag(u)**2)) * dV 
#endif 
                end do
!$OMP END PARALLEL DO
#endif 
#elif(_DIM_==3)
                n1 = 1;  if (m%g%n1min==0) n1 = 2
                n2 = 1;  if (m%g%n2min==0) n2 = 2
                n3 = 1;  if (m%g%n3min==0) n3 = 2
                nv1 = lbound(m%V,1) + n1 - 1
                nv2 = lbound(m%V,2) + n2 - 1
                nv3 = lbound(m%V,3) + n3 - 1
#ifndef _OPENMP
#ifdef _REAL_
                E_pot = sum(m%V(nv1:,nv2:,nv3:) * this%u(n1:,n2:,n3:)**2) * dV 
#else
                E_pot = sum(m%V(nv1:,nv2:,nv3:) * (real(this%u(n1:,n2:,n3:), prec)**2 &
                                                 +aimag(this%u(n1:,n2:,n3:))**2)) * dV 
#endif 
#else
                E_pot = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, V) REDUCTION(+:E_pot)
                do j=1,n_threads
                      u => this%u(n1:,n2:,max(n3,lbound(this%u,3)+m%g%jj(j-1)):lbound(this%u,3)+m%g%jj(j)-1)
                      V => m%V(nv1:,nv2:,max(nv3,lbound(m%V,3)+m%g%jj(j-1)):lbound(m%V,3)+m%g%jj(j)-1)
#ifdef _REAL_
                      E_pot = E_pot + sum(V * u**2) * dV 
#else
                      E_pot = E_pot + sum(V * (real(u, prec)**2+aimag(u)**2)) * dV 
#endif 
                end do
!$OMP END PARALLEL DO
#endif 
#endif
            end select      
#endif
      
#ifdef _MPI_
            call MPI_Reduce(E_pot, h, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            E_pot = h
            call MPI_Bcast(E_pot, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif            
        else
            E_pot = 0.0_prec
        end if    

        class default
           stop "E: wrong spectral method for schroedinger wave function"
        end select
    end function potential_energy

    
    function interaction_energy(this) result(E_int)
        class(S(wf_schroedinger)), intent(inout) :: this
        real(kind=prec) :: E_int
        real(kind=prec) :: dV 
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: h 
#endif
        integer :: n1, n2, n3
#ifdef _OPENMP
#if(_DIM_==1)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:)
        real(kind=prec), pointer :: V(:)
#elif(_DIM_==2)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:)
        real(kind=prec), pointer :: V(:,:)
#elif(_DIM_==3)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:,:)
        real(kind=prec), pointer :: V(:,:,:)
#endif
#if defined(_HERMITE_)
        real(kind=prec), pointer :: w(:)
#endif        
        integer :: j
#endif        


        select type (m=>this%m); class is (S(schroedinger))

        if (m%cubic_coupling/=0_prec) then
            call this%to_real_space

#if defined(_HERMITE_)
#ifndef _OPENMP
#ifdef _REAL_
            E_int = 0.5_prec*m%cubic_coupling*sum(this%u**4 &
#else
            E_int = 0.5_prec*m%cubic_coupling*sum((real(this%u, prec)**2+imag(this%u)**2)**2 &
#endif
#if(_DIM_==1)
                *m%g%weights_x) 
#elif(_DIM_==2)
                *spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1) &
                *spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1))
#elif(_DIM_==3)
                *spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), 3, m%g%n3max-m%g%n3min+1) &
                *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), 3, m%g%n3max-m%g%n3min+1) &
                *spread(spread(m%g%weights_z, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) )
#endif
#else
          E_int = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, w) REDUCTION(+:E_int)
          do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                w => m%g%weights_x(lbound(m%g%weights_x,1)+m%g%jj(j-1):lbound(m%g%weights_x,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                w => m%g%weights_y(lbound(m%g%weights_y,1)+m%g%jj(j-1):lbound(m%g%weights_y,1)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                w => m%g%weights_z(lbound(m%g%weights_z,1)+m%g%jj(j-1):lbound(m%g%weights_z,1)+m%g%jj(j)-1)
#endif 
#ifdef _REAL_
               E_int = E_int + 0.5_prec*m%cubic_coupling*sum(u**4 &
#else
               E_int = E_int + 0.5_prec*m%cubic_coupling*sum((real(u, prec)**2+imag(u)**2)**2 &
#endif
#if(_DIM_==1)
                *w) 
#elif(_DIM_==2)
                *spread(m%g%weights_x, 2, m%g%jj(j)-m%g%jj(j-1)) &
                *spread(w, 1, m%g%n1max-m%g%n1min+1))
#elif(_DIM_==3)
                *spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), 3, m%g%jj(j)-m%g%jj(j-1)) &
                *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), 3, m%g%jj(j)-m%g%jj(j-1)) &
                *spread(spread(w, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) )
#endif

          end do
!$OMP END PARALLEL DO
#endif 

#else

#if(_DIM_==1)
            dV = m%g%dx 
#elif(_DIM_==2)
            dV = m%g%dx * m%g%dy 
#elif(_DIM_==3)
            dV = m%g%dx * m%g%dy * m%g%dz
#endif

            select case(m%boundary_conditions)
            case(periodic, dirichlet)
#ifndef _OPENMP
#ifdef _REAL_
                E_int = 0.5_prec*m%cubic_coupling*sum(this%u**4) * dV
#else
                E_int = 0.5_prec*m%cubic_coupling*sum((real(this%u, prec)**2+imag(this%u)**2)**2) * dV
#endif
#else
          E_int = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u) REDUCTION(+:E_int)
          do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
#endif
#ifdef _REAL_
                E_int = E_int + 0.5_prec*m%cubic_coupling*sum(u**4) * dV
#else
                E_int = E_int + 0.5_prec*m%cubic_coupling*sum((real(u, prec)**2+imag(u)**2)**2) * dV
#endif
          end do
!$OMP END PARALLEL DO
#endif 
            case(neumann) ! for summation exlude first index in each dimension 
#if(_DIM_==1)
                n1 = 1;  if (m%g%n1min==0) n1 = 2
#ifndef _OPENMP
#ifdef _REAL_
                E_int = 0.5_prec*m%cubic_coupling*sum(this%u(n1:)**4) * dV
#else
                E_int = 0.5_prec*m%cubic_coupling*sum((real(this%u(n1:), prec)**2+aimag(this%u(n1:)**2))**2) * dV
#endif
#else
                E_int = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u) REDUCTION(+:E_int)
                do j=1,n_threads
                      u => this%u(max(n1,lbound(this%u,1)+m%g%jj(j-1)):lbound(this%u,1)+m%g%jj(j)-1)
#ifdef _REAL_
                      E_int = E_int + 0.5_prec*m%cubic_coupling*sum(u**4) * dV
#else
                      E_int = E_int + 0.5_prec*m%cubic_coupling*sum((real(u, prec)**2+aimag(u)**2)**2) * dV
#endif 
                end do
!$OMP END PARALLEL DO
#endif 
#elif(_DIM_==2)
                n1 = 1;  if (m%g%n1min==0) n1 = 2
                n2 = 1;  if (m%g%n2min==0) n2 = 2
#ifndef _OPENMP
#ifdef _REAL_
                E_int = 0.5_prec*m%cubic_coupling*sum(this%u(n1:,n2:)**4) * dV
#else
                E_int = 0.5_prec*m%cubic_coupling*sum((real(this%u(n1:,n2:), prec)**2 &
                                                      +aimag(this%u(n1:,n2:))**2)**2) * dV
#endif
#else
                E_int = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u) REDUCTION(+:E_int)
                do j=1,n_threads
                      u => this%u(n1:,max(n2,lbound(this%u,2)+m%g%jj(j-1)):lbound(this%u,2)+m%g%jj(j)-1)
#ifdef _REAL_
                      E_int = E_int + 0.5_prec*m%cubic_coupling*sum(u**4) * dV
#else
                      E_int = E_int + 0.5_prec*m%cubic_coupling*sum((real(u, prec)**2+aimag(u)**2)**2) * dV
#endif 
                end do
!$OMP END PARALLEL DO
#endif 
#elif(_DIM_==3)
                n1 = 1;  if (m%g%n1min==0) n1 = 2
                n2 = 1;  if (m%g%n2min==0) n2 = 2
                n3 = 1;  if (m%g%n3min==0) n3 = 2
#ifndef _OPENMP
#ifdef _REAL_
                E_int = 0.5_prec*m%cubic_coupling*sum(this%u(n1:,n2:,n3:)**4) * dV
#else
                E_int = 0.5_prec*m%cubic_coupling*sum((real(this%u(n1:,n2:,n3:), prec)**2 &
                                                       +aimag(this%u(n1:,n2:,n3:))**2)**2) * dV
#endif
#else
                E_int = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u) REDUCTION(+:E_int)
                do j=1,n_threads
                      u => this%u(n1:,n2:,max(n3,lbound(this%u,3)+m%g%jj(j-1)):lbound(this%u,3)+m%g%jj(j)-1)
#ifdef _REAL_
                      E_int = E_int + 0.5_prec*m%cubic_coupling*sum(u**4) * dV
#else
                      E_int = E_int + 0.5_prec*m%cubic_coupling*sum((real(u, prec)**2+aimag(u)**2)**2) * dV
#endif
                end do
!$OMP END PARALLEL DO
#endif 
#endif
            end select       

#endif     

#ifdef _MPI_
            call MPI_Reduce(E_int, h, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            E_int = h
            call MPI_Bcast(E_int, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif            
        else
            E_int = 0.0_prec
        end if    

        class default
           stop "E: wrong spectral method for schroedinger wave function"
        end select
    end function interaction_energy



#if(_DIM_==1)
   subroutine get_realspace_observables(this, E_pot, E_int, x_mean, x_dev)
#elif(_DIM_==2)
   subroutine get_realspace_observables(this, E_pot, E_int, x_mean, x_dev, y_mean, y_dev)
#elif(_DIM_==3)
   subroutine get_realspace_observables(this, E_pot, E_int, x_mean, x_dev, y_mean, y_dev, z_mean, z_dev)
#endif
        class(S(wf_schroedinger)), intent(inout) :: this
        real(kind=prec), intent(out) :: E_pot
        real(kind=prec), intent(out) :: E_int
        real(kind=prec), intent(out) :: x_mean, x_dev
#if(_DIM_>=2)
        real(kind=prec), intent(out) :: y_mean, y_dev
#endif        
#if(_DIM_>=3)
        real(kind=prec), intent(out) :: z_mean, z_dev
#endif        
        real(kind=prec) :: dV 
        integer :: n1, n2, n3
        integer :: nv1, nv2, nv3
#ifdef _OPENMP
#if(_DIM_==1)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:)
        real(kind=prec), pointer :: V(:)
#elif(_DIM_==2)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:)
        real(kind=prec), pointer :: V(:,:)
#elif(_DIM_==3)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: u(:,:,:)
        real(kind=prec), pointer :: V(:,:,:)
#endif
#if defined(_HERMITE_)
        real(kind=prec), pointer :: w(:)
#endif        
        real(kind=prec), pointer :: nodes(:)
        integer :: j
#endif   
#ifdef _MPI_        
        integer :: ierr
        real(kind=prec) :: buf(2+_DIM_) 
        real(kind=prec) :: buf1(2+_DIM_) 
#endif        



        select type (m=>this%m ); class is (S(schroedinger))

        call this%to_real_space

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(_HERMITE_)
#ifndef _OPENMP
#if(_DIM_==1)
          x_mean = sum(m%g%weights_x*m%g%nodes_x &
#elif(_DIM_==2)
          x_mean = sum(spread(m%g%weights_x*m%g%nodes_x, 2, m%g%n2max-m%g%n2min+1) &
                      *spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1) &
#elif(_DIM_==3)
          x_mean = sum(spread(spread(m%g%weights_x*m%g%nodes_x, 2, m%g%n2max-m%g%n2min+1), &
                   3, m%g%n3max-m%g%n3min+1) &
                     *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), 3, m%g%n3max-m%g%n3min+1) &
                     *spread(spread(m%g%weights_z, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) &
#endif
#ifdef _REAL_
                     *(this%u**2) )
#else
                     *(real(this%u,kind=prec)**2 +aimag(this%u)**2) )
#endif 
#else
          x_mean = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, w, nodes) REDUCTION(+:x_mean)
          do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                w => m%g%weights_x(lbound(m%g%weights_x,1)+m%g%jj(j-1):lbound(m%g%weights_x,1)+m%g%jj(j)-1)
                nodes => m%g%nodes_x(lbound(m%g%nodes_x,1)+m%g%jj(j-1):lbound(m%g%nodes_x,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                w => m%g%weights_y(lbound(m%g%weights_y,1)+m%g%jj(j-1):lbound(m%g%weights_y,1)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                w => m%g%weights_z(lbound(m%g%weights_z,1)+m%g%jj(j-1):lbound(m%g%weights_z,1)+m%g%jj(j)-1)
#endif 
#if(_DIM_==1)
                x_mean = x_mean + sum(w*nodes &
#elif(_DIM_==2)
                x_mean = x_mean + sum(spread(m%g%weights_x*m%g%nodes_x, 2, m%g%jj(j)-m%g%jj(j-1)) &
                      *spread(w, 1, m%g%n1max-m%g%n1min+1) &
#elif(_DIM_==3)
                x_mean = x_mean + sum(spread(spread(m%g%weights_x*m%g%nodes_x, 2, m%g%n2max-m%g%n2min+1), &
                   3, m%g%jj(j)-m%g%jj(j-1)) &
                     *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), 3, m%g%jj(j)-m%g%jj(j-1) ) &
                     *spread(spread(w, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) &
#endif
#ifdef _REAL_
                     *(u**2) )
#else
                     *(real(u,kind=prec)**2 +aimag(u)**2) )
#endif 

          end do
!$OMP END PARALLEL DO
#endif 
#if(_DIM_>=2)
#ifndef _OPENMP
#if(_DIM_==2)
          y_mean = sum(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1) &
                      *spread(m%g%weights_y*m%g%nodes_y, 1, m%g%n1max-m%g%n1min+1) &
#elif(_DIM_==3)
          y_mean = sum(spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), &
                   3, m%g%n3max-m%g%n3min+1) &
                     *spread(spread(m%g%weights_y*m%g%nodes_y, 1, m%g%n1max-m%g%n1min+1), &
                   3, m%g%n3max-m%g%n3min+1) &
                     *spread(spread(m%g%weights_z, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) &
#endif
#ifdef _REAL_
                     *(this%u**2) )
#else
                     *(real(this%u,kind=prec)**2 +aimag(this%u)**2) )
#endif 
#else
          y_mean = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, w, nodes) REDUCTION(+:y_mean)
          do j=1,n_threads
#if(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                w => m%g%weights_y(lbound(m%g%weights_y,1)+m%g%jj(j-1):lbound(m%g%weights_y,1)+m%g%jj(j)-1)
                nodes => m%g%nodes_y(lbound(m%g%nodes_y,1)+m%g%jj(j-1):lbound(m%g%nodes_y,1)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                w => m%g%weights_z(lbound(m%g%weights_z,1)+m%g%jj(j-1):lbound(m%g%weights_z,1)+m%g%jj(j)-1)
#endif 
#if(_DIM_==2)
                y_mean = y_mean + sum(spread(m%g%weights_x, 2, m%g%jj(j)-m%g%jj(j-1)) &
                      *spread(w*nodes, 1, m%g%n1max-m%g%n1min+1) &
#elif(_DIM_==3)
                y_mean =  y_mean + sum(spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), &
                   3, m%g%jj(j)-m%g%jj(j-1) ) &
                     *spread(spread(m%g%weights_y*m%g%nodes_y, 1, m%g%n1max-m%g%n1min+1), &
                   3, m%g%jj(j)-m%g%jj(j-1) ) &
                     *spread(spread(w, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) &
#endif
#ifdef _REAL_
                     *(u**2) )
#else
                     *(real(u,kind=prec)**2 +aimag(u)**2) )
#endif 

          end do
!$OMP END PARALLEL DO
#endif 
#endif 

#if(_DIM_==3)
#ifndef _OPENMP
          z_mean = sum(spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), &
                   3, m%g%n3max-m%g%n3min+1) &
                     *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), &
                   3, m%g%n3max-m%g%n3min+1) &
                     *spread(spread(m%g%weights_z*m%g%nodes_z, 1, m%g%n1max-m%g%n1min+1), &
                   2, m%g%n2max-m%g%n2min+1) &
#ifdef _REAL_
                     *(this%u**2) )
#else
                     *(real(this%u,kind=prec)**2 +aimag(this%u)**2) )
#endif 
#else
          z_mean = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, w, nodes) REDUCTION(+:z_mean)
          do j=1,n_threads
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                w => m%g%weights_z(lbound(m%g%weights_z,1)+m%g%jj(j-1):lbound(m%g%weights_z,1)+m%g%jj(j)-1)
                nodes => m%g%nodes_z(lbound(m%g%nodes_z,1)+m%g%jj(j-1):lbound(m%g%nodes_z,1)+m%g%jj(j)-1)
                z_mean = z_mean + sum(spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), &
                   3, m%g%jj(j)-m%g%jj(j-1)) &
                     *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), &
                   3, m%g%jj(j)-m%g%jj(j-1)) &
                     *spread(spread(w*nodes, 1, m%g%n1max-m%g%n1min+1), &
                   2, m%g%n2max-m%g%n2min+1) &
#ifdef _REAL_
                     *(u**2) )
#else
                     *(real(u,kind=prec)**2 +aimag(u)**2) )
#endif 
          end do
!$OMP END PARALLEL DO
#endif 
#endif

#else

#if(_DIM_==1)
        dV = m%g%dx 
#elif(_DIM_==2)
        dV = m%g%dx * m%g%dy 
#elif(_DIM_==3)
        dV = m%g%dx * m%g%dy * m%g%dz
#endif
        select case(m%boundary_conditions)
        case(periodic, dirichlet)
#ifndef _OPENMP
#if(_DIM_==1)
            x_mean = sum( m%g%nodes_x & 
#elif(_DIM_==2)
            x_mean = sum(spread(m%g%nodes_x, 2, m%g%n2max-m%g%n2min+1) &
#elif(_DIM_==3)
            x_mean =sum(spread(spread(m%g%nodes_x, 2, m%g%n2max-m%g%n2min+1), &
                3, m%g%n3max-m%g%n3min+1) &
#endif
#ifdef _REAL_
                * this%u**2) * dV 
#else
                * (real(this%u, prec)**2+aimag(this%u)**2)) * dV 
#endif 
#else
            x_mean = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, nodes) REDUCTION(+:x_mean)
            do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                nodes => m%g%nodes_x(lbound(m%g%nodes_x,1)+m%g%jj(j-1):lbound(m%g%nodes_x,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
#endif 
#if(_DIM_==1)
                x_mean = x_mean + sum( nodes & 
#elif(_DIM_==2)
                x_mean = x_mean + sum(spread(m%g%nodes_x, 2, m%g%jj(j)-m%g%jj(j-1) ) &
#elif(_DIM_==3)
                x_mean =x_mean + sum(spread(spread(m%g%nodes_x, 2, m%g%n2max-m%g%n2min+1), &
                3, m%g%jj(j)-m%g%jj(j-1) ) &
#endif
#ifdef _REAL_
                * u**2) * dV 
#else
                * (real(u, prec)**2+aimag(u)**2)) * dV 
#endif 
            end do
!$OMP END PARALLEL DO
#endif 

#if(_DIM_>=2)
#ifndef _OPENMP
#if(_DIM_==2)
            y_mean = sum(spread(m%g%nodes_y, 1, m%g%n1max-m%g%n1min+1) &
#elif(_DIM_==3)
            y_mean =sum(spread(spread(m%g%nodes_y, 1, m%g%n1max-m%g%n1min+1), &
                3, m%g%n3max-m%g%n3min+1) &
#endif
#ifdef _REAL_
                * this%u**2) * dV 
#else
                * (real(this%u, prec)**2+aimag(this%u)**2)) * dV 
#endif 
#else
            y_mean = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, nodes) REDUCTION(+:y_mean)
            do j=1,n_threads
#if(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                nodes => m%g%nodes_y(lbound(m%g%nodes_y,1)+m%g%jj(j-1):lbound(m%g%nodes_y,1)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
#endif 
#if(_DIM_==2)
                y_mean = y_mean + sum(spread(nodes, 1, m%g%n1max-m%g%n1min+1) &
#elif(_DIM_==3)
                y_mean = y_mean + sum(spread(spread(m%g%nodes_y, 1, m%g%n1max-m%g%n1min+1), &
                                3, m%g%jj(j)-m%g%jj(j-1)) &
#endif
#ifdef _REAL_
                * u**2) * dV 
#else
                * (real(u, prec)**2+aimag(u)**2)) * dV 
#endif 

            end do
!$OMP END PARALLEL DO
#endif 
#endif 

#if(_DIM_==3)
#ifndef _OPENMP
            z_mean = sum(spread(spread(m%g%nodes_z, 1, m%g%n1max-m%g%n1min+1), &
                2, m%g%n2max-m%g%n2min+1) &
#ifdef _REAL_
                * this%u**2) * dV 
#else
                * (real(this%u, prec)**2+aimag(this%u)**2)) * dV 
#endif 
#else
             z_mean = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, nodes) REDUCTION(+:z_mean)
             do j=1,n_threads
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                nodes => m%g%nodes_z(lbound(m%g%nodes_z,1)+m%g%jj(j-1):lbound(m%g%nodes_z,1)+m%g%jj(j)-1)
                z_mean = z_mean + sum(spread(spread(nodes, 1, m%g%n1max-m%g%n1min+1), &
                                  2, m%g%n2max-m%g%n2min+1) &
#ifdef _REAL_
                * u**2) * dV 
#else
                * (real(u, prec)**2+aimag(u)**2)) * dV 
#endif 

             end do
!$OMP END PARALLEL DO
#endif 
#endif 

        case(neumann) ! for summation exlude first index in each dimension 
!!!TODO: x2,y2,z2 for Neumann !!!!

#if(_DIM_==1)
            n1 = 1;  if (m%g%n1min==0) n1 = 2
#elif(_DIM_==2)
            n1 = 1;  if (m%g%n1min==0) n1 = 2
            n2 = 1;  if (m%g%n2min==0) n2 = 2
#elif(_DIM_==3)
            n1 = 1;  if (m%g%n1min==0) n1 = 2
            n2 = 1;  if (m%g%n2min==0) n2 = 2
            n3 = 1;  if (m%g%n3min==0) n3 = 2
#endif
         end select      

#endif
      
#ifdef _MPI_
            buf(1) = x_mean
#if(_DIM_>=2)
            buf(2) = y_mean
#endif            
#if(_DIM_>=3)
            buf(3) = z_mean
#endif            
            call MPI_Reduce(buf, buf1, _DIM_, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(buf1, _DIM_, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            x_mean = buf1(1)
#if(_DIM_>=2)
            y_mean = buf1(2)
#endif            
#if(_DIM_>=3)
            z_mean = buf1(3)
#endif            
#endif            



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        E_pot = 0.0_prec
        E_int = 0.0_prec


#if defined(_HERMITE_)
#ifndef _OPENMP
#if(_DIM_==1)
          x_dev = sum(m%g%weights_x*(m%g%nodes_x-x_mean)**2 &
#elif(_DIM_==2)
          x_dev = sum(spread(m%g%weights_x*(m%g%nodes_x-x_mean)**2, 2, m%g%n2max-m%g%n2min+1) &
                      *spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1) &
#elif(_DIM_==3)
          x_dev = sum(spread(spread(m%g%weights_x*(m%g%nodes_x-x_mean)**2, 2, m%g%n2max-m%g%n2min+1), &
                   3, m%g%n3max-m%g%n3min+1) &
                     *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), 3, m%g%n3max-m%g%n3min+1) &
                     *spread(spread(m%g%weights_z, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) &
#endif
#ifdef _REAL_
                     *(this%u**2) )
#else
                     *(real(this%u,kind=prec)**2 +aimag(this%u)**2) )
#endif 
#else
          x_dev = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, w, nodes) REDUCTION(+:x_dev)
          do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                w => m%g%weights_x(lbound(m%g%weights_x,1)+m%g%jj(j-1):lbound(m%g%weights_x,1)+m%g%jj(j)-1)
                nodes => m%g%nodes_x(lbound(m%g%nodes_x,1)+m%g%jj(j-1):lbound(m%g%nodes_x,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                w => m%g%weights_y(lbound(m%g%weights_y,1)+m%g%jj(j-1):lbound(m%g%weights_y,1)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                w => m%g%weights_z(lbound(m%g%weights_z,1)+m%g%jj(j-1):lbound(m%g%weights_z,1)+m%g%jj(j)-1)
#endif 
#if(_DIM_==1)
                x_dev = x_dev + sum(w*(nodes-x_mean)**2 &
#elif(_DIM_==2)
                x_dev = x_dev + sum(spread(m%g%weights_x*(m%g%nodes_x-x_mean)**2, 2, m%g%jj(j)-m%g%jj(j-1)) &
                      *spread(w, 1, m%g%n1max-m%g%n1min+1) &
#elif(_DIM_==3)
               x_dev = x_dev + sum(spread(spread(m%g%weights_x*(m%g%nodes_x-x_mean)**2, 2, m%g%n2max-m%g%n2min+1), &
                   3, m%g%jj(j)-m%g%jj(j-1)) &
                     *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), 3, m%g%jj(j)-m%g%jj(j-1)) &
                     *spread(spread(w, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) &
#endif
#ifdef _REAL_
                     *(u**2) )
#else
                     *(real(u,kind=prec)**2 +aimag(u)**2) )
#endif 

          end do
!$OMP END PARALLEL DO
#endif 
#if(_DIM_>=2)
#ifndef _OPENMP
#if(_DIM_==2)
          y_dev = sum(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1) &
                      *spread(m%g%weights_y*(m%g%nodes_y-y_mean)**2, 1, m%g%n1max-m%g%n1min+1) &
#elif(_DIM_==3)
          y_dev = sum(spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), &
                   3, m%g%n3max-m%g%n3min+1) &
                     *spread(spread(m%g%weights_y*(m%g%nodes_y-y_mean)**2, 1, m%g%n1max-m%g%n1min+1), &
                   3, m%g%n3max-m%g%n3min+1) &
                     *spread(spread(m%g%weights_z, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) &
#endif
#ifdef _REAL_
                     *(this%u**2) )
#else
                     *(real(this%u,kind=prec)**2 +aimag(this%u)**2) )
#endif 
#else
          y_dev = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, w, nodes) REDUCTION(+:y_dev)
          do j=1,n_threads
#if(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                w => m%g%weights_y(lbound(m%g%weights_y,1)+m%g%jj(j-1):lbound(m%g%weights_y,1)+m%g%jj(j)-1)
                nodes => m%g%nodes_y(lbound(m%g%nodes_y,1)+m%g%jj(j-1):lbound(m%g%nodes_y,1)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                w => m%g%weights_z(lbound(m%g%weights_z,1)+m%g%jj(j-1):lbound(m%g%weights_z,1)+m%g%jj(j)-1)
#endif 
#if(_DIM_==2)
                y_dev = y_dev + sum(spread(m%g%weights_x, 2, m%g%jj(j)-m%g%jj(j-1)) &
                      *spread(w*(nodes-y_mean)**2, 1, m%g%n1max-m%g%n1min+1) &
#elif(_DIM_==3)
                y_dev = y_dev + sum(spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), &
                   3, m%g%jj(j)-m%g%jj(j-1)) &
                     *spread(spread(m%g%weights_y*(m%g%nodes_y-y_mean)**2, 1, m%g%n1max-m%g%n1min+1), &
                   3, m%g%jj(j)-m%g%jj(j-1)) &
                     *spread(spread(w, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) &
#endif
#ifdef _REAL_
                     *(u**2) )
#else
                     *(real(u,kind=prec)**2 +aimag(u)**2) )
#endif 
          end do
!$OMP END PARALLEL DO
#endif 
#endif 

#if(_DIM_==3)
#ifndef _OPENMP
          z_dev = sum(spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), &
                   3, m%g%n3max-m%g%n3min+1) &
                     *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), &
                   3, m%g%n3max-m%g%n3min+1) &
                     *spread(spread(m%g%weights_z*(m%g%nodes_z-z_mean)**2, 1, m%g%n1max-m%g%n1min+1), &
                   2, m%g%n2max-m%g%n2min+1) &
#ifdef _REAL_
                     *(this%u**2) )
#else
                     *(real(this%u,kind=prec)**2 +aimag(this%u)**2) )
#endif 
#else
          z_dev = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, w, nodes) REDUCTION(+:z_dev)
          do j=1,n_threads
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                w => m%g%weights_z(lbound(m%g%weights_z,1)+m%g%jj(j-1):lbound(m%g%weights_z,1)+m%g%jj(j)-1)
                nodes => m%g%nodes_z(lbound(m%g%nodes_z,1)+m%g%jj(j-1):lbound(m%g%nodes_z,1)+m%g%jj(j)-1)
                z_dev = z_dev + sum(spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), &
                   3, m%g%jj(j)-m%g%jj(j-1)) &
                     *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), &
                   3, m%g%jj(j)-m%g%jj(j-1)) &
                     *spread(spread(w*(nodes-z_mean)**2, 1, m%g%n1max-m%g%n1min+1), &
                   2, m%g%n2max-m%g%n2min+1) &
#ifdef _REAL_
                     *(u**2) )
#else
                     *(real(u,kind=prec)**2 +aimag(u)**2) )
#endif 

          end do
!$OMP END PARALLEL DO
#endif 
#endif

        if (associated(m%V)) then
#ifndef _OPENMP
#if(_DIM_==1)
          E_pot = sum(m%g%weights_x &
#elif(_DIM_==2)
          E_pot = sum(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1) &
                      *spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1) &
#elif(_DIM_==3)
          E_pot = sum(spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), 3, m%g%n3max-m%g%n3min+1) &
                     *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), 3, m%g%n3max-m%g%n3min+1) &
                     *spread(spread(m%g%weights_z, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) &
#endif
#ifdef _REAL_
                     *m%V *(this%u**2) )
#else
                     *m%V *(real(this%u,kind=prec)**2 +aimag(this%u)**2) )
#endif 
#else
          E_pot = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, V, w) REDUCTION(+:E_pot)
          do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                V => m%V(lbound(m%V,1)+m%g%jj(j-1):lbound(m%V,1)+m%g%jj(j)-1)
                w => m%g%weights_x(lbound(m%g%weights_x,1)+m%g%jj(j-1):lbound(m%g%weights_x,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                V => m%V(:,lbound(m%V,2)+m%g%jj(j-1):lbound(m%V,2)+m%g%jj(j)-1)
                w => m%g%weights_y(lbound(m%g%weights_y,1)+m%g%jj(j-1):lbound(m%g%weights_y,1)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                V => m%V(:,:,lbound(m%V,3)+m%g%jj(j-1):lbound(m%V,3)+m%g%jj(j)-1)
                w => m%g%weights_z(lbound(m%g%weights_z,1)+m%g%jj(j-1):lbound(m%g%weights_z,1)+m%g%jj(j)-1)
#endif 
#if(_DIM_==1)
                E_pot = E_pot + sum(w &
#elif(_DIM_==2)
                E_pot = E_pot + sum(spread(m%g%weights_x, 2, m%g%jj(j)-m%g%jj(j-1)) &
                      *spread(w, 1, m%g%n1max-m%g%n1min+1) &
#elif(_DIM_==3)
                E_pot = E_pot + sum(spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), 3, m%g%jj(j)-m%g%jj(j-1)) &
                     *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), 3, m%g%jj(j)-m%g%jj(j-1)) &
                     *spread(spread(w, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) &
#endif
#ifdef _REAL_
                     *V *(u**2) )
#else
                     *V *(real(u,kind=prec)**2 +aimag(u)**2) )
#endif 
          end do
!$OMP END PARALLEL DO
#endif 
        endif

        if (m%cubic_coupling/=0_prec) then
#ifndef _OPENMP
#ifdef _REAL_
               E_int = 0.5_prec*m%cubic_coupling*sum(this%u**4 &
#else
               E_int = 0.5_prec*m%cubic_coupling*sum((real(this%u, prec)**2+imag(this%u)**2)**2 &
#endif
#if(_DIM_==1)
                *m%g%weights_x) 
#elif(_DIM_==2)
                *spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1) &
                *spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1))
#elif(_DIM_==3)
                *spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), 3, m%g%n3max-m%g%n3min+1) &
                *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), 3, m%g%n3max-m%g%n3min+1) &
                *spread(spread(m%g%weights_z, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) )
#endif
#else
          E_int = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, w) REDUCTION(+:E_int)
          do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                w => m%g%weights_x(lbound(m%g%weights_x,1)+m%g%jj(j-1):lbound(m%g%weights_x,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                w => m%g%weights_y(lbound(m%g%weights_y,1)+m%g%jj(j-1):lbound(m%g%weights_y,1)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                w => m%g%weights_z(lbound(m%g%weights_z,1)+m%g%jj(j-1):lbound(m%g%weights_z,1)+m%g%jj(j)-1)
#endif 
#ifdef _REAL_
               E_int = E_int + 0.5_prec*m%cubic_coupling*sum(u**4 &
#else
               E_int = E_int + 0.5_prec*m%cubic_coupling*sum((real(u, prec)**2+imag(u)**2)**2 &
#endif
#if(_DIM_==1)
                *w) 
#elif(_DIM_==2)
                *spread(m%g%weights_x, 2, m%g%jj(j)-m%g%jj(j-1)) &
                *spread(w, 1, m%g%n1max-m%g%n1min+1))
#elif(_DIM_==3)
                *spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), 3, m%g%jj(j)-m%g%jj(j-1)) &
                *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), 3, m%g%jj(j)-m%g%jj(j-1)) &
                *spread(spread(w, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) )
#endif

          end do
!$OMP END PARALLEL DO
#endif 
        endif
        
        
#else

        select case(m%boundary_conditions)
        case(periodic, dirichlet)
#ifndef _OPENMP
#if(_DIM_==1)
            x_dev = sum( (m%g%nodes_x-x_mean)**2 & 
#elif(_DIM_==2)
            x_dev = sum(spread((m%g%nodes_x-x_mean)**2, 2, m%g%n2max-m%g%n2min+1) &
#elif(_DIM_==3)
            x_dev =sum(spread(spread((m%g%nodes_x-x_mean)**2, 2, m%g%n2max-m%g%n2min+1), &
                3, m%g%n3max-m%g%n3min+1) &
#endif
#ifdef _REAL_
                * this%u**2) * dV 
#else
                * (real(this%u, prec)**2+aimag(this%u)**2)) * dV 
#endif 
#else
            x_dev = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, nodes) REDUCTION(+:x_dev)
            do j=1,n_threads
#if(_DIM_==1)
                u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                nodes => m%g%nodes_x(lbound(m%g%nodes_x,1)+m%g%jj(j-1):lbound(m%g%nodes_x,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
#endif 
#if(_DIM_==1)
                x_dev = x_dev + sum( (nodes-x_mean)**2 & 
#elif(_DIM_==2)
                x_dev = x_dev + sum(spread((m%g%nodes_x-x_mean)**2, 2, m%g%jj(j)-m%g%jj(j-1)) &
#elif(_DIM_==3)
                x_dev =x_dev + sum(spread(spread((m%g%nodes_x-x_mean)**2, 2, m%g%n2max-m%g%n2min+1), &
                3, m%g%jj(j)-m%g%jj(j-1)) &
#endif
#ifdef _REAL_
                * u**2) * dV 
#else
                * (real(u, prec)**2+aimag(u)**2)) * dV 
#endif 
            end do
!$OMP END PARALLEL DO
#endif 

#if(_DIM_>=2)
#ifndef _OPENMP
#if(_DIM_==2)
            y_dev = sum(spread((m%g%nodes_y-y_mean)**2, 1, m%g%n1max-m%g%n1min+1) &
#elif(_DIM_==3)
            y_dev =sum(spread(spread((m%g%nodes_y-y_mean)**2, 1, m%g%n1max-m%g%n1min+1), &
                3, m%g%n3max-m%g%n3min+1) &
#endif
#ifdef _REAL_
                * this%u**2) * dV 
#else
                * (real(this%u, prec)**2+aimag(this%u)**2)) * dV 
#endif 
#else
            y_dev = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, nodes) REDUCTION(+:y_dev)
            do j=1,n_threads
#if(_DIM_==2)
                u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                nodes => m%g%nodes_y(lbound(m%g%nodes_y,1)+m%g%jj(j-1):lbound(m%g%nodes_y,1)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
#endif 
#if(_DIM_==2)
                y_dev = y_dev + sum(spread((nodes-y_mean)**2, 1, m%g%n1max-m%g%n1min+1) &
#elif(_DIM_==3)
                y_dev = y_dev + sum(spread(spread((m%g%nodes_y-y_mean)**2, 1, m%g%n1max-m%g%n1min+1), &
                                    3, m%g%jj(j)-m%g%jj(j-1) ) &
#endif
#ifdef _REAL_
                * u**2) * dV 
#else
                * (real(u, prec)**2+aimag(u)**2)) * dV 
#endif 

            end do
!$OMP END PARALLEL DO
#endif 
#endif 

#if(_DIM_==3)
#ifndef _OPENMP
            z_dev =sum(spread(spread((m%g%nodes_z-z_mean)**2, 1, m%g%n1max-m%g%n1min+1), &
                2, m%g%n2max-m%g%n2min+1) &
#ifdef _REAL_
                * this%u**2) * dV 
#else
                * (real(this%u, prec)**2+aimag(this%u)**2)) * dV 
#endif 
#else
             z_dev = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, nodes) REDUCTION(+:z_dev)
             do j=1,n_threads
                u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                nodes => m%g%nodes_z(lbound(m%g%nodes_z,1)+m%g%jj(j-1):lbound(m%g%nodes_z,1)+m%g%jj(j)-1)
                z_dev = z_dev + sum(spread(spread((nodes-z_mean)**2, 1, m%g%n1max-m%g%n1min+1), &
                   2, m%g%n2max-m%g%n2min+1) &
#ifdef _REAL_
                   * u**2) * dV 
#else
                   * (real(u, prec)**2+aimag(u)**2)) * dV 
#endif 
             end do
!$OMP END PARALLEL DO
#endif 
#endif 

            if (associated(m%V)) then
#ifndef _OPENMP
#ifdef _REAL_
                E_pot = sum(m%V * this%u**2) * dV 
#else
                E_pot = sum(m%V * (real(this%u, prec)**2+aimag(this%u)**2)) * dV 
#endif 
#else
                E_pot = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, V) REDUCTION(+:E_pot)
                do j=1,n_threads
#if(_DIM_==1)
                      u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                      V => m%V(lbound(m%V,1)+m%g%jj(j-1):lbound(m%V,1)+m%g%jj(j)-1)
      #elif(_DIM_==2)
                      u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                      V => m%V(:,lbound(m%V,2)+m%g%jj(j-1):lbound(m%V,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                      u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                      V => m%V(:,:,lbound(m%V,3)+m%g%jj(j-1):lbound(m%V,3)+m%g%jj(j)-1)
#endif
#ifdef _REAL_
                      E_pot = E_pot + sum(V * u**2) * dV 
#else
                      E_pot = E_pot + sum(V * (real(u, prec)**2+aimag(u)**2)) * dV 
#endif 
                end do
!$OMP END PARALLEL DO
#endif 
            endif  
            if (m%cubic_coupling/=0_prec) then
#ifndef _OPENMP
#ifdef _REAL_
                E_int = 0.5_prec*m%cubic_coupling*sum(this%u**4) * dV
#else
                E_int = 0.5_prec*m%cubic_coupling*sum((real(this%u, prec)**2+imag(this%u)**2)**2) * dV
#endif
#else
                E_int = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u) REDUCTION(+:E_int)
                do j=1,n_threads
#if(_DIM_==1)
                   u => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                   u => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                   u => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
#endif
#ifdef _REAL_
                   E_int = E_int + 0.5_prec*m%cubic_coupling*sum(u**4) * dV
#else
                   E_int = E_int + 0.5_prec*m%cubic_coupling*sum((real(u, prec)**2+imag(u)**2)**2) * dV
#endif
                end do
!$OMP END PARALLEL DO
#endif 
            endif
        case(neumann) ! for summation exlude first index in each dimension 
!!!TODO: x2,y2,z2 for Neumann !!!!

#if(_DIM_==1)
            if (associated(m%V)) then
                nv1 = lbound(m%V,1) + n1 - 1
#ifndef _OPENMP
#ifdef _REAL_
                E_pot = sum(m%V(nv1:) * this%u(n1:)**2) * dV 
#else
                E_pot = sum(m%V(nv1:) * (real(this%u(n1:), prec)**2+aimag(this%u(n1:))**2)) * dV 
#endif 
#else
                E_pot = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, V) REDUCTION(+:E_pot)
                do j=1,n_threads
                      u => this%u(max(n1,lbound(this%u,1)+m%g%jj(j-1)):lbound(this%u,1)+m%g%jj(j)-1)
                      V => m%V(max(nv1,lbound(m%V,1)+m%g%jj(j-1)):lbound(m%V,1)+m%g%jj(j)-1)
#ifdef _REAL_
                      E_pot = E_pot + sum(V * u**2) * dV 
#else
                      E_pot = E_pot + sum(V * (real(u, prec)**2+aimag(u)**2)) * dV 
#endif 
                end do
!$OMP END PARALLEL DO
#endif 
            endif  
            if (m%cubic_coupling/=0_prec) then
#ifndef _OPENMP
#ifdef _REAL_
                E_int = 0.5_prec*m%cubic_coupling*sum(this%u(n1:)**4) * dV
#else
                E_int = 0.5_prec*m%cubic_coupling*sum((real(this%u(n1:), prec)**2+aimag(this%u(n1:)**2))**2) * dV
#endif
#else
                E_int = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u) REDUCTION(+:E_int)
                do j=1,n_threads
                      u => this%u(max(n1,lbound(this%u,1)+m%g%jj(j-1)):lbound(this%u,1)+m%g%jj(j)-1)
#ifdef _REAL_
                      E_int = E_int + 0.5_prec*m%cubic_coupling*sum(u**4) * dV
#else
                      E_int = E_int + 0.5_prec*m%cubic_coupling*sum((real(u, prec)**2+aimag(u)**2)**2) * dV
#endif 
                end do
!$OMP END PARALLEL DO
#endif 
            endif
#elif(_DIM_==2)
            if (associated(m%V)) then
                nv1 = lbound(m%V,1) + n1 - 1
                nv2 = lbound(m%V,2) + n2 - 1
#ifndef _OPENMP
#ifdef _REAL_
                E_pot = sum(m%V(nv1:,nv2:) * this%u(n1:,n2:)**2) * dV 
#else
                E_pot = sum(m%V(nv1:,nv2:) * (real(this%u(n1:,n2:), prec)**2 &
                                            +aimag(this%u(n1:,n2:))**2)) * dV 
#endif 
#else
                E_pot = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, V) REDUCTION(+:E_pot)
                do j=1,n_threads
                      u => this%u(n1:,max(n2,lbound(this%u,2)+m%g%jj(j-1)):lbound(this%u,2)+m%g%jj(j)-1)
                      V => m%V(nv1:,max(nv2,lbound(m%V,2)+m%g%jj(j-1)):lbound(m%V,2)+m%g%jj(j)-1)
#ifdef _REAL_
                      E_pot = E_pot + sum(V * u**2) * dV 
#else
                      E_pot = E_pot + sum(V * (real(u, prec)**2+aimag(u)**2)) * dV 
#endif 
                end do
!$OMP END PARALLEL DO
#endif 
            endif  
            if (m%cubic_coupling/=0_prec) then
#ifndef _OPENMP
#ifdef _REAL_
                E_int = 0.5_prec*m%cubic_coupling*sum(this%u(n1:,n2:)**4) * dV
#else
                E_int = 0.5_prec*m%cubic_coupling*sum((real(this%u(n1:,n2:), prec)**2 &
                                                      +aimag(this%u(n1:,n2:))**2)**2) * dV
#endif
#else
                E_int = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u) REDUCTION(+:E_int)
                do j=1,n_threads
                      u => this%u(n1:,max(n2,lbound(this%u,2)+m%g%jj(j-1)):lbound(this%u,2)+m%g%jj(j)-1)
#ifdef _REAL_
                      E_int = E_int + 0.5_prec*m%cubic_coupling*sum(u**4) * dV
#else
                      E_int = E_int + 0.5_prec*m%cubic_coupling*sum((real(u, prec)**2+aimag(u)**2)**2) * dV
#endif 
                end do
!$OMP END PARALLEL DO
#endif 
            endif
#elif(_DIM_==3)
            if (associated(m%V)) then
                nv1 = lbound(m%V,1) + n1 - 1
                nv2 = lbound(m%V,2) + n2 - 1
                nv3 = lbound(m%V,3) + n3 - 1
#ifndef _OPENMP
#ifdef _REAL_
                E_pot = sum(m%V(nv1:,nv2:,nv3:) * this%u(n1:,n2:,n3:)**2) * dV 
#else
                E_pot = sum(m%V(nv1:,nv2:,nv3:) * (real(this%u(n1:,n2:,n3:), prec)**2 &
                                                 +aimag(this%u(n1:,n2:,n3:))**2)) * dV 
#endif 
#else
                E_pot = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u, V) REDUCTION(+:E_pot)
                do j=1,n_threads
                      u => this%u(n1:,n2:,max(n3,lbound(this%u,3)+m%g%jj(j-1)):lbound(this%u,3)+m%g%jj(j)-1)
                      V => m%V(nv1:,nv2:,max(nv3,lbound(m%V,3)+m%g%jj(j-1)):lbound(m%V,3)+m%g%jj(j)-1)
#ifdef _REAL_
                      E_pot = E_pot + sum(V * u**2) * dV 
#else
                      E_pot = E_pot + sum(V * (real(u, prec)**2+aimag(u)**2)) * dV 
#endif 
                end do
!$OMP END PARALLEL DO
#endif 
            endif  
            if (m%cubic_coupling/=0_prec) then
#ifndef _OPENMP
#ifdef _REAL_
                E_int = 0.5_prec*m%cubic_coupling*sum(this%u(n1:,n2:,n3:)**4) * dV
#else
                E_int = 0.5_prec*m%cubic_coupling*sum((real(this%u(n1:,n2:,n3:), prec)**2 &
                                                       +aimag(this%u(n1:,n2:,n3:))**2)**2) * dV
#endif
#else
                E_int = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u) REDUCTION(+:E_int)
                do j=1,n_threads
                      u => this%u(n1:,n2:,max(n3,lbound(this%u,3)+m%g%jj(j-1)):lbound(this%u,3)+m%g%jj(j)-1)
#ifdef _REAL_
                      E_int = E_int + 0.5_prec*m%cubic_coupling*sum(u**4) * dV
#else
                      E_int = E_int + 0.5_prec*m%cubic_coupling*sum((real(u, prec)**2+aimag(u)**2)**2) * dV
#endif
                end do
!$OMP END PARALLEL DO
#endif 
            endif
#endif
            end select      

#endif
      
#ifdef _MPI_
            buf(1) = E_pot
            buf(2) = E_int
            buf(3) = x_dev
#if(_DIM_>=2)
            buf(4) = y_dev
#endif            
#if(_DIM_>=3)
            buf(5) = z_dev
#endif            
            call MPI_Reduce(buf, buf1, 2+_DIM_, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            call MPI_Bcast(buf1, 2+_DIM_, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
            E_pot = buf1(1)
            E_int = buf1(2)
            x_dev = buf1(3)
#if(_DIM_>=2)
            y_dev = buf1(4)
#endif            
#if(_DIM_>=3)
            z_dev = buf1(5)
#endif            
#endif            
            x_dev = sqrt(x_dev)
#if(_DIM_>=2)
            y_dev = sqrt(y_dev)
#endif            
#if(_DIM_>=3)
            z_dev = sqrt(z_dev)
#endif            
        class default
           stop "E: wrong spectral method for schroedinger wave function"
        end select
   end subroutine get_realspace_observables



    function observable(this, f) result(obs)
        class(S(wf_schroedinger)), intent(inout) :: this
        real(kind=prec), external :: f
        real(kind=prec) :: obs
#if(_DIM_==1)
        real(kind=prec) :: x
        integer :: ix
#elif(_DIM_==2)
        real(kind=prec) :: x, y
        integer :: ix, iy
#elif(_DIM_==3)
        real(kind=prec) :: x, y, z
        integer :: ix, iy, iz
#endif
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: h 
#endif
        integer :: n1, n2, n3
        integer :: ox, oy, oz

        select type (m=>this%m ); class is (S(schroedinger))

        call this%to_real_space

        n1 = m%g%n1min
        ox = lbound(this%u, 1) - n1 
#if(_DIM_>=2)
        n2 = m%g%n2min
        oy = lbound(this%u, 2) - n2 
#endif
#if(_DIM_>=3)
        n3 = m%g%n3min
        oz = lbound(this%u, 3) - n3 
#endif

#if defined(_HERMITE_)

        obs = 0.0_prec
#if(_DIM_==1)
!$OMP PARALLEL DO PRIVATE(ix) REDUCTION(+:obs)
        do ix = n1+ox,m%g%n1max+ox
            obs = obs + f(m%g%nodes_x(ix-ox))*m%g%weights_x(ix-ox) &
#ifdef _REAL_
                   *this%u(ix)**2
#else
                   *(real(this%u(ix),prec)**2+aimag(this%u(ix))**2)
#endif
        end do
!$OMP END PARALLEL DO 
#elif(_DIM_==2)       
!$OMP PARALLEL DO PRIVATE(ix, iy) REDUCTION(+:obs)
        do iy = n2+oy, m%g%n2max+oy
            do ix = n1+ox,m%g%n1max+ox
                obs = obs + f(m%g%nodes_x(ix-ox),m%g%nodes_x(iy-oy)) &
                      *m%g%weights_x(ix-ox)*m%g%weights_y(iy-oy) &
#ifdef _REAL_
                      *this%u(ix,iy)**2
#else
                      *(real(this%u(ix,iy),prec)**2+aimag(this%u(ix,iy))**2)
#endif
            end do
        end do     
!$OMP END PARALLEL DO 
#elif(_DIM_==3)
!$OMP PARALLEL DO PRIVATE(ix, iy, iz) REDUCTION(+:obs)
        do iz = n3+oz, m%g%n3max+oz
           do iy = n2+oy, m%g%n2max+oy
                do ix = n1+ox,m%g%n1max+ox
                    obs = obs + f(m%g%nodes_x(ix-ox),m%g%nodes_y(iy-oy),m%g%nodes_z(iz-oz)) &
                          *m%g%weights_x(ix-ox)*m%g%weights_y(iy-oy)*m%g%weights_z(iz-oz) &
#ifdef _REAL_
                          *this%u(ix,iy,iz)**2
#else
                          *(real(this%u(ix,iy,iz),prec)**2+aimag(this%u(ix,iy,iz))**2)
#endif
                end do
           end do     
        end do    
!$OMP END PARALLEL DO 
#endif

#else
        select case(this%m%boundary_conditions)
        case(periodic, dirichlet)
        case(neumann) ! for summation exlude first index in each dimension 
            if (n1==0) n1 = 1 
#if(_DIM_>=2)
            if (n2==0) n2 = 1
#endif
#if(_DIM_>=3)
            if (n3==0) n3 = 1
#endif
        end select

        obs = 0.0_prec
#if(_DIM_==1)
!$OMP PARALLEL DO PRIVATE(ix, x) REDUCTION(+:obs)
        do ix = n1+ox,m%g%n1max+ox
            x = m%g%xmin + m%g%dx*(ix-ox)    
#ifdef _REAL_
            obs = obs + f(x)*this%u(ix)**2
#else
            obs = obs + f(x)*(real(this%u(ix),prec)**2+aimag(this%u(ix))**2)
#endif
        end do
!$OMP END PARALLEL DO 
        obs = obs * m%g%dx 
#elif(_DIM_==2)       
!$OMP PARALLEL DO PRIVATE(ix, x, iy, y) REDUCTION(+:obs)
        do iy = n2+oy, m%g%n2max+oy
            y = m%g%xmin + m%g%dy*(iy-oy)    
            do ix = n1+ox,m%g%n1max+ox
                x = m%g%xmin + m%g%dx*(ix-ox)    
#ifdef _REAL_
                obs = obs + f(x,y)*this%u(ix,iy)**2
#else
                obs = obs + f(x,y)*(real(this%u(ix,iy),prec)**2+aimag(this%u(ix,iy))**2)
#endif
            end do
        end do     
!$OMP END PARALLEL DO 
        obs = obs * m%g%dx * m%g%dy
#elif(_DIM_==3)
!$OMP PARALLEL DO PRIVATE(ix, x, iy, y, iz, z) REDUCTION(+:obs)
        do iz = n3+oz, m%g%n3max+oz
           z = m%g%zmin + m%g%dz*(iz-oz)    
           do iy = n2+oy, m%g%n2max+oy
                y = m%g%xmin + m%g%dy*(iy-oy)    
                do ix = n1+ox,m%g%n1max+ox
                    x = m%g%xmin + m%g%dx*(ix-ox)    
#ifdef _REAL_
                    obs = obs + f(x,y,z)*this%u(ix,iy,iz)**2
#else
                    obs = obs + f(x,y,z)*(real(this%u(ix,iy,iz),prec)**2+aimag(this%u(ix,iy,iz))**2)
#endif
                end do
           end do     
        end do    
!$OMP END PARALLEL DO 
        obs = obs * m%g%dx * m%g%dy * m%g%dz 
#endif

#endif 

#ifdef _MPI_
            call MPI_Reduce(obs, h, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            obs = h
            call MPI_Bcast(obs, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif            
        class default
           stop "E: wrong spectral method for schroedinger wave function"
        end select
    end function observable



    subroutine get_energy_expectation_deviation(this, E_expectation, E_deviation)
       class(S(wf_schroedinger)), intent(inout) :: this
       real(kind=prec), intent(out) :: E_expectation
       real(kind=prec), intent(out) :: E_deviation
#ifdef _MPI_
       integer :: ierr
       real(kind=prec) :: h
#endif
       !real(kind=prec) :: E2
       real(kind=prec) :: E_var, dV
       integer :: n1, n2, n3
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
#if defined(_HERMITE_)
        real(kind=prec), pointer :: w(:)
#endif        
        integer :: j
#endif        

       select type (m=>this%m ); class is (S(schroedinger))

       call m%initialize_tmp

       m%tmp%u = 0.0_prec
       call this%add_apply_A(m%tmp)
       call this%add_apply_B(m%tmp)
       call this%to_real_space
       call m%tmp%to_real_space

#if defined(_HERMITE_)

#ifndef _OPENMP
#ifdef _REAL_
       E_expectation = -sum(m%tmp%u * this%u &
#else
       E_expectation = sum(aimag(conjg(m%tmp%u)*this%u) &
#endif
#if(_DIM_==1)
                *m%g%weights_x) 
#elif(_DIM_==2)
                *spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1) &
                *spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1))
#elif(_DIM_==3)
                *spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), 3, m%g%n3max-m%g%n3min+1) &
                *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), 3, m%g%n3max-m%g%n3min+1) &
                *spread(spread(m%g%weights_z, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) )
#endif
#else
       E_expectation = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u1, u2, w) REDUCTION(+:E_expectation)
            do j=1,n_threads
#if(_DIM_==1)
                u1 => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                u2 => m%tmp%u(lbound(m%tmp%u,1)+m%g%jj(j-1):lbound(m%tmp%u,1)+m%g%jj(j)-1)
                w => m%g%weights_x(lbound(m%g%weights_x,1)+m%g%jj(j-1):lbound(m%g%weights_x,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u1 => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                u2 => m%tmp%u(:,lbound(m%tmp%u,2)+m%g%jj(j-1):lbound(m%tmp%u,2)+m%g%jj(j)-1)
                w => m%g%weights_y(lbound(m%g%weights_y,1)+m%g%jj(j-1):lbound(m%g%weights_y,1)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u1 => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                u2 => m%tmp%u(:,:,lbound(m%tmp%u,3)+m%g%jj(j-1):lbound(m%tmp%u,3)+m%g%jj(j)-1)
                w => m%g%weights_z(lbound(m%g%weights_z,1)+m%g%jj(j-1):lbound(m%g%weights_z,1)+m%g%jj(j)-1)
#endif 
#ifdef _REAL_
                E_expectation = E_expectation + (-1.0_prec)*sum(u2 * u1 &
#else
                E_expectation = E_expectation + sum(aimag(conjg(u2)* u1) &
#endif
#if(_DIM_==1)
                   *w) 
#elif(_DIM_==2)
                   *spread(m%g%weights_x, 2, m%g%jj(j)-m%g%jj(j-1) ) &
                   *spread(w, 1, m%g%n1max-m%g%n1min+1))
#elif(_DIM_==3)
                   *spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), 3, m%g%jj(j)-m%g%jj(j-1) ) &
                   *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), 3, m%g%jj(j)-m%g%jj(j-1) ) &
                   *spread(spread(w, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) )
#endif

            end do
!$OMP END PARALLEL DO
#endif         

#else

#if(_DIM_==1)
       dV = m%g%dx
#elif(_DIM_==2)
       dV = m%g%dx*m%g%dy
#elif(_DIM_==3)
       dV = m%g%dx*m%g%dy*m%g%dz
#endif

       select case(m%boundary_conditions)
       case(periodic, dirichlet)
           ! TODO check negative signs
#ifndef _OPENMP
#ifdef _REAL_
           E_expectation = -sum(m%tmp%u * this%u)*dV
#else
           E_expectation = sum(aimag(conjg(m%tmp%u)*this%u))*dV
#endif 
#else
           E_expectation = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u1, u2) REDUCTION(+:E_expectation)
            do j=1,n_threads
#if(_DIM_==1)
                u1 => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                u2 => m%tmp%u(lbound(m%tmp%u,1)+m%g%jj(j-1):lbound(m%tmp%u,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u1 => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                u2 => m%tmp%u(:,lbound(m%tmp%u,2)+m%g%jj(j-1):lbound(m%tmp%u,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u1 => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                u2 => m%tmp%u(:,:,lbound(m%tmp%u,3)+m%g%jj(j-1):lbound(m%tmp%u,3)+m%g%jj(j)-1)
#endif 
#ifdef _REAL_
                E_expectation = E_expectation + (-1.0_prec)* sum(u2 * u1)*dV
#else
                E_expectation = E_expectation + sum(aimag(conjg(u2)*u1))*dV
#endif 
            end do
!$OMP END PARALLEL DO
#endif     
       case(neumann) ! for summation exlude first index in each dimension 
#if(_DIM_==1)
           n1 = 1;  if (m%g%n1min==0) n1 = 2
#ifndef _OPENMP
#ifdef _REAL_
           E_expectation = -sum(m%tmp%u(n1:) * this%u(n1:))*dV
#else
           E_expectation = sum(aimag(conjg(m%tmp%u(n1:))*this%u(n1:)))*dV
#endif 
#else
           E_expectation = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u1, u2) REDUCTION(+:E_expectation)
           do j=1,n_threads
               u1 => this%u(max(n1,lbound(this%u,1)+m%g%jj(j-1)):lbound(this%u,1)+m%g%jj(j)-1)
               u2 => m%tmp%u(max(n1,lbound(m%tmp%u,1)+m%g%jj(j-1)):lbound(m%tmp%u,1)+m%g%jj(j)-1)
#ifdef _REAL_
               E_expectation = E_expectation + (-1.0_prec)* sum(u2 * u1)*dV
#else
               E_expectation = E_expectation + sum(aimag(conjg(u2)*u1))*dV
#endif 
           end do
!$OMP END PARALLEL DO
#endif 
#elif(_DIM_==2)
           n1 = 1;  if (m%g%n1min==0) n1 = 2
           n2 = 1;  if (m%g%n2min==0) n2 = 2
#ifndef _OPENMP
#ifdef _REAL_
           E_expectation = -sum(m%tmp%u(n1:,n2:) * this%u(n1:,n2:))*dV
#else
           E_expectation = sum(aimag(conjg(m%tmp%u(n1:,n2:))*this%u(n1:,n2:)))*dV
#endif 
#else
           E_expectation = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u1, u2) REDUCTION(+:E_expectation)
           do j=1,n_threads
               u1 => this%u(n1:,max(n2,lbound(this%u,2)+m%g%jj(j-1)):lbound(this%u,2)+m%g%jj(j)-1)
               u2 => m%tmp%u(n1:,max(n2,lbound(m%tmp%u,2)+m%g%jj(j-1)):lbound(m%tmp%u,2)+m%g%jj(j)-1)
#ifdef _REAL_
               E_expectation = E_expectation + (-1.0_prec)* sum(u2 * u1)*dV
#else
               E_expectation = E_expectation + sum(aimag(conjg(u2)*u1))*dV
#endif 
           end do
!$OMP END PARALLEL DO
#endif 
#elif(_DIM_==3)
           n1 = 1;  if (m%g%n1min==0) n1 = 2
           n2 = 1;  if (m%g%n2min==0) n2 = 2
           n3 = 1;  if (m%g%n3min==0) n3 = 2
#ifndef _OPENMP
#ifdef _REAL_
           E_expectation = -sum(m%tmp%u(n1:,n2:,n3:) * this%u(n1:,n2:,n3:))*dV
#else
           E_expectation = sum(aimag(conjg(m%tmp%u(n1:,n2:,n3:))*this%u(n1:,n2:,n3:)))*dV
#endif 
#else
           E_expectation = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u1, u2) REDUCTION(+:E_expectation)
           do j=1,n_threads
               u1 => this%u(n1:,n2:,max(n3,lbound(this%u,3)+m%g%jj(j-1)):lbound(this%u,3)+m%g%jj(j)-1)
               u2 => m%tmp%u(n1:,n2:,max(n3,lbound(m%tmp%u,3)+m%g%jj(j-1)):lbound(m%tmp%u,3)+m%g%jj(j)-1)
#ifdef _REAL_
               E_expectation = E_expectation + (-1.0_prec)* sum(u2 * u1)*dV
#else
               E_expectation = E_expectation + sum(aimag(conjg(u2)*u1))*dV
#endif 
           end do
!$OMP END PARALLEL DO
#endif 
#endif
       end select

#endif

#ifdef _MPI_
       call MPI_Reduce(E_expectation, h, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(h, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       E_expectation = h 
#endif

#if defined(_HERMITE_)

#ifndef _OPENMP
#ifdef _REAL_
       E_var = sum( (-m%tmp%u - E_expectation*this%u )**2  &
#else
       E_var = sum( ((real(m%tmp%u, prec) - E_expectation*aimag(this%u) )**2 &
                +   (aimag(m%tmp%u) + E_expectation*real(this%u, prec))**2) &
#endif
#if(_DIM_==1)
                *m%g%weights_x) 
#elif(_DIM_==2)
                *spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1) &
                *spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1))
#elif(_DIM_==3)
                *spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), 3, m%g%n3max-m%g%n3min+1) &
                *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), 3, m%g%n3max-m%g%n3min+1) &
                *spread(spread(m%g%weights_z, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) )
#endif
#else
       E_var = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u1, u2, w) REDUCTION(+:E_var)
            do j=1,n_threads
#if(_DIM_==1)
                u1 => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                u2 => m%tmp%u(lbound(m%tmp%u,1)+m%g%jj(j-1):lbound(m%tmp%u,1)+m%g%jj(j)-1)
                w => m%g%weights_x(lbound(m%g%weights_x,1)+m%g%jj(j-1):lbound(m%g%weights_x,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u1 => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                u2 => m%tmp%u(:,lbound(m%tmp%u,2)+m%g%jj(j-1):lbound(m%tmp%u,2)+m%g%jj(j)-1)
                w => m%g%weights_y(lbound(m%g%weights_y,1)+m%g%jj(j-1):lbound(m%g%weights_y,1)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u1 => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                u2 => m%tmp%u(:,:,lbound(m%tmp%u,3)+m%g%jj(j-1):lbound(m%tmp%u,3)+m%g%jj(j)-1)
                w => m%g%weights_z(lbound(m%g%weights_z,1)+m%g%jj(j-1):lbound(m%g%weights_z,1)+m%g%jj(j)-1)
#endif 
#ifdef _REAL_
                E_var = E_var + sum( (-u2 - E_expectation*u1 )**2  &
#else
                E_var = E_var + sum( ((real(u2, prec) - E_expectation*aimag(u1) )**2 &
                +   (aimag(u2) + E_expectation*real(u1, prec))**2) &
#endif
#if(_DIM_==1)
                *w) 
#elif(_DIM_==2)
                *spread(m%g%weights_x, 2, m%g%jj(j)-m%g%jj(j-1)) &
                *spread(w, 1, m%g%n1max-m%g%n1min+1))
#elif(_DIM_==3)
                *spread(spread(m%g%weights_x, 2, m%g%n2max-m%g%n2min+1), 3, m%g%jj(j)-m%g%jj(j-1)) &
                *spread(spread(m%g%weights_y, 1, m%g%n1max-m%g%n1min+1), 3, m%g%jj(j)-m%g%jj(j-1)) &
                *spread(spread(w, 1, m%g%n1max-m%g%n1min+1), 2, m%g%n2max-m%g%n2min+1) )
#endif
            end do
!$OMP END PARALLEL DO
#endif  

#else

       select case(m%boundary_conditions)
       case(periodic, dirichlet)
           ! TODO check negative signs
#ifndef _OPENMP
#ifdef _REAL_
           E_var = sum( (-m%tmp%u - E_expectation*this%u )**2 )*dV
#else
           E_var = sum( (real(m%tmp%u, prec) - E_expectation*aimag(this%u) )**2 &
                    +   (aimag(m%tmp%u) + E_expectation*real(this%u, prec))**2)*dV
#endif
#else
           E_var = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u1, u2) REDUCTION(+:E_var)
            do j=1,n_threads
#if(_DIM_==1)
                u1 => this%u(lbound(this%u,1)+m%g%jj(j-1):lbound(this%u,1)+m%g%jj(j)-1)
                u2 => m%tmp%u(lbound(m%tmp%u,1)+m%g%jj(j-1):lbound(m%tmp%u,1)+m%g%jj(j)-1)
#elif(_DIM_==2)
                u1 => this%u(:,lbound(this%u,2)+m%g%jj(j-1):lbound(this%u,2)+m%g%jj(j)-1)
                u2 => m%tmp%u(:,lbound(m%tmp%u,2)+m%g%jj(j-1):lbound(m%tmp%u,2)+m%g%jj(j)-1)
#elif(_DIM_==3)
                u1 => this%u(:,:,lbound(this%u,3)+m%g%jj(j-1):lbound(this%u,3)+m%g%jj(j)-1)
                u2 => m%tmp%u(:,:,lbound(m%tmp%u,3)+m%g%jj(j-1):lbound(m%tmp%u,3)+m%g%jj(j)-1)
#endif 
#ifdef _REAL_
                E_var = E_var + sum( (-u2 - E_expectation*u1 )**2 )*dV
#else
                E_var = E_var + sum( (real(u2, prec) - E_expectation*aimag(u1) )**2 &
                         +   (aimag(u2) + E_expectation*real(u1, prec))**2)*dV
#endif
            end do
!$OMP END PARALLEL DO
#endif  
       case(neumann) ! for summation exlude first index in each dimension 
#if(_DIM_==1)
           n1 = 1;  if (m%g%n1min==0) n1 = 2
#ifndef _OPENMP
#ifdef _REAL_
           E_var = sum( (-m%tmp%u(n1:) - E_expectation*this%u(n1:) )**2 )*dV
#else
           E_var = sum( (real(m%tmp%u(n1:), prec) - E_expectation*real(this%u(n1:), prec) )**2 &
                    +   (aimag(m%tmp%u(n1:))-E_expectation*aimag(this%u(n1:)))**2)*dV
#endif
#else
           E_var = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u1, u2) REDUCTION(+:E_var)
           do j=1,n_threads
               u1 => this%u(max(n1,lbound(this%u,1)+m%g%jj(j-1)):lbound(this%u,1)+m%g%jj(j)-1)
               u2 => m%tmp%u(max(n1,lbound(m%tmp%u,1)+m%g%jj(j-1)):lbound(m%tmp%u,1)+m%g%jj(j)-1)
#ifdef _REAL_
               E_var = E_var + sum( (-u2 - E_expectation*u1 )**2 )*dV
#else
               E_var = E_var + sum( (real(u2, prec) - E_expectation*aimag(u1) )**2 &
                         +   (aimag(u2) + E_expectation*real(u1, prec))**2)*dV
#endif 
           end do
!$OMP END PARALLEL DO
#endif 
#elif(_DIM_==2)
           n1 = 1;  if (m%g%n1min==0) n1 = 2
           n2 = 1;  if (m%g%n2min==0) n2 = 2
#ifndef _OPENMP
#ifdef _REAL_
           E_var = sum( (-m%tmp%u(n1:,n2:) - E_expectation*this%u(n1:,n2:) )**2 )*dV
#else
           E_var = sum( (real(m%tmp%u(n1:,n2:), prec) - E_expectation*real(this%u(n1:,n2:), prec) )**2 &
                    +   (aimag(m%tmp%u(n1:,n2:))-E_expectation*aimag(this%u(n1:,n2:)))**2)*dV
#endif
#else
           E_var = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u1, u2) REDUCTION(+:E_var)
           do j=1,n_threads
               u1 => this%u(n1:,max(n2,lbound(this%u,2)+m%g%jj(j-1)):lbound(this%u,2)+m%g%jj(j)-1)
               u2 => m%tmp%u(n1:,max(n2,lbound(m%tmp%u,2)+m%g%jj(j-1)):lbound(m%tmp%u,2)+m%g%jj(j)-1)
#ifdef _REAL_
               E_var = E_var + sum( (-u2 - E_expectation*u1 )**2 )*dV
#else
               E_var = E_var + sum( (real(u2, prec) - E_expectation*aimag(u1) )**2 &
                         +   (aimag(u2) + E_expectation*real(u1, prec))**2)*dV
#endif 
           end do
!$OMP END PARALLEL DO
#endif 
#elif(_DIM_==3)
           n1 = 1;  if (m%g%n1min==0) n1 = 2
           n2 = 1;  if (m%g%n2min==0) n2 = 2
           n3 = 1;  if (m%g%n3min==0) n3 = 2
#ifndef _OPENMP
#ifdef _REAL_
           E_var = sum( (-m%tmp%u(n1:,n2:,n3:) - E_expectation*this%u(n1:,n2:,n3:) )**2 )*dV
#else
           E_var = sum( (real(m%tmp%u(n1:,n2:,n3:), prec) - E_expectation*real(this%u(n1:,n2:,n3:), prec) )**2 &
                    +   (aimag(m%tmp%u(n1:,n2:,n3:))-E_expectation*aimag(this%u(n1:,n2:,n3:)))**2)*dV
#endif
#else
           E_var = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, u1, u2) REDUCTION(+:E_var)
           do j=1,n_threads
               u1 => this%u(n1:,n2:,max(n3,lbound(this%u,3)+m%g%jj(j-1)):lbound(this%u,3)+m%g%jj(j)-1)
               u2 => m%tmp%u(n1:,n2:,max(n3,lbound(m%tmp%u,3)+m%g%jj(j-1)):lbound(m%tmp%u,3)+m%g%jj(j)-1)
#ifdef _REAL_
               E_var = E_var + sum( (-u2 - E_expectation*u1 )**2 )*dV
#else
               E_var = E_var + sum( (real(u2, prec) - E_expectation*aimag(u1) )**2 &
                         +   (aimag(u2) + E_expectation*real(u1, prec))**2)*dV
#endif 
           end do
!$OMP END PARALLEL DO
#endif 
#endif
       end select

#endif

#ifdef _MPI_
       call MPI_Reduce(E_var, h, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(h, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       E_var = h 
#endif
       E_deviation = m%hbar*sqrt(E_var)
       E_expectation = m%hbar*E_expectation 

       class default
           stop "E: wrong spectral method for schroedinger wave function"
       end select
    end subroutine get_energy_expectation_deviation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    subroutine solution_out(psi, t, dt, flag, abort_flag)
        class(_WAVE_FUNCTION_), intent(inout) :: psi
        real(kind=prec), intent(in) :: t 
        real(kind=prec), intent(in) :: dt 
        character, intent(in) :: flag
        logical, intent(out) :: abort_flag
        integer, save :: step
        real(kind=prec) :: norm, E_kin, E_pot, E_int, E_tot, x_mean, x_dev
#if (_DIM_>=2)
        real(kind=prec) :: y_mean, y_dev, r_dev 
#endif 
#if (_DIM_>=3)
        real(kind=prec) :: z_mean, z_dev
#endif 
        integer(kind=8), save :: calc 
        real :: calc_time         

        select type (psi); class is (S(wf_schroedinger))

        if (flag=="i" .or. flag=="I") then ! init =======================
            call tick(calc)
            step = 0
            if (this_proc==0) then
                write (*,*), "#   #                   t                  dt                norm" & 
            //"               E_kin               E_pot               E_int               E_tot" & 
#if(_DIM_==1)
            //"              x_mean               x_dev cumm.time"
#elif(_DIM_==2)
           // "              x_mean               x_dev              y_mean               y_dev" & 
           // "               r_dev   cummulated time" 
#elif(_DIM_==3)
           // "              x_mean               x_dev              y_mean               y_dev" & 
           // "              z_mean               z_dev               r_dev   cummulated time" 
#endif 
                write (*,*), "#---------------------------------------------------------------------------" &
                       //"------------------------------------------------------------------" &
                       //"------------------------------------------------------" 
            end if           
        end if   

        if (flag=="i" .or. flag=="I" .or. flag=="e" .or. flag=="E") then ! execute =====
            norm = psi%norm()
            E_kin = psi%kinetic_energy()
#if(_DIM_==1)
            call psi%get_realspace_observables(E_pot, E_int, x_mean, x_dev)
#elif(_DIM_==2)
            call psi%get_realspace_observables(E_pot, E_int, x_mean, x_dev, y_mean, y_dev)
            r_dev = sqrt(x_dev**2 + y_dev**2)
#elif(_DIM_==3)
            call psi%get_realspace_observables(E_pot, E_int, x_mean, x_dev, y_mean, y_dev, z_mean, z_dev)
            r_dev = sqrt(x_dev**2 + y_dev**2 + z_dev**2)
#endif
            E_tot = E_kin + E_pot + E_int
            calc_time = tock(calc)
            if (this_proc==0) then
#if(_DIM_==1)
              write (*, '(I6,9E20.12,F10.2)') step, t,  dt, norm, E_kin, E_pot, E_int, E_tot, x_mean, x_dev, calc_time 
#elif(_DIM_==2)
              write (*, '(I6,12E20.12,F10.2)') step, t,  dt, norm, E_kin, E_pot, E_int, E_tot, x_mean, x_dev, &
                                             y_mean, y_dev, r_dev, calc_time 
#elif(_DIM_==3)
              write (*, '(I6,14E20.12,F10.2)') step, t,  dt, norm, E_kin, E_pot, E_int, E_tot, x_mean, x_dev, &
                                             y_mean, y_dev, z_mean, z_dev, r_dev, calc_time 
#endif
            end if
            step = step + 1
        end if
    
        if (flag=="f" .or. flag=="F") then ! finish ====================
        end if   

        end select

        abort_flag = .false.
    end subroutine solution_out


    subroutine extrapolation_imaginary_time_step(this, dt, extrapolation_order, psi0, psi1, start_with_B)
        class(S(wf_schroedinger)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt
        integer, intent(in) :: extrapolation_order
        class(S(wf_schroedinger)), intent(inout) :: psi0, psi1
        logical, intent(in), optional :: start_with_B

        logical :: start_with_B_1
        integer :: k, j 
       _COMPLEX_OR_REAL_(kind=prec) :: dt1
        real(kind=prec) :: c(extrapolation_order)

        start_with_B_1 = .false.
        if (present(start_with_B)) then
            start_with_B_1 = start_with_B
        end if

        select case(extrapolation_order)
        case (1)
           c(1) = 1.0_prec
        case (2)
           c(1) = -1.0_prec/3.0_prec
           c(2) =  4.0_prec/3.0_prec
        case (3)
           c(1) =  1.0_prec/24.0_prec
           c(2) = -16.0_prec/15.0_prec
           c(3) =  81.0_prec/40.0_prec
        case (4)
           c(1) = -1.0_prec/360.0_prec
           c(2) =  16.0_prec/45.0_prec
           c(3) = -729.0_prec/280.0_prec
           c(4) =  1024.0_prec/315.0_prec
        case (5)
           c(1) =  1.0_prec/8640
           c(2) =  -64.0_prec/945.0_prec
           c(3) =  6561.0_prec/4480.0_prec
           c(4) = -16384.0_prec/2835.0_prec
           c(5) =  390625.0_prec/72576.0_prec
        case (6)
           c(1) = -1.0_prec/302400
           c(2) = 8.0_prec/945.0_prec
           c(3) = -2187.0_prec/4480.0_prec
           c(4) = 65536.0_prec/14175.0_prec
           c(5) = -9765625.0_prec/798336.0_prec
           c(6) = 17496.0_prec/1925.0_prec
        case default
           print *, "Error: extrapolation_order must be in 1..6"
           stop
        end select

        if (start_with_B_1) then
           call this%to_real_space()
           call psi0%copy(this)
           call psi1%copy(this)

           call this%imaginary_time_propagate_B(0.5_prec*dt, 0)
           call this%imaginary_time_propagate_A(dt)
           call this%imaginary_time_propagate_B(0.5_prec*dt, 1)

           this%u = c(1)*this%u

           do k = 2, extrapolation_order
               dt1 = dt/real(k, kind=prec)
               call psi1%copy(psi0)
               call psi1%imaginary_time_propagate_B(0.5_prec*dt1, 0)
               do j = 1, k-1
                   call psi1%imaginary_time_propagate_A(dt1)
                   call psi1%imaginary_time_propagate_B(dt1, 2)
               end do    
               call psi1%imaginary_time_propagate_A(dt1)
               call psi1%imaginary_time_propagate_B(0.5_prec*dt1, 1)

               this%u = this%u + c(k)*psi1%u
           end do
        else
           call this%to_frequency_space()
           call psi0%copy(this)
           call psi1%copy(this)

           call this%imaginary_time_propagate_A(0.5_prec*dt)
           call this%imaginary_time_propagate_B(dt, 2)
           call this%imaginary_time_propagate_A(0.5_prec*dt)
           
           this%u = c(1)*this%u

           do k = 2, extrapolation_order
               dt1 = dt/real(k, kind=prec)
               call psi1%copy(psi0)
               call psi1%imaginary_time_propagate_A(0.5_prec*dt1)
               do j = 1, k-1
                   call psi1%imaginary_time_propagate_B(dt1, 2)
                   call psi1%imaginary_time_propagate_A(dt1)
               end do
               call psi1%imaginary_time_propagate_B(dt1, 2)
               call psi1%imaginary_time_propagate_A(0.5_prec*dt1)

               this%u = this%u + c(k)*psi1%u
           end do
        end if

    end subroutine extrapolation_imaginary_time_step
    



    subroutine splitting_imaginary_time_step(this, dt, splitting_scheme, start_with_B)
        class(S(wf_schroedinger)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: splitting_scheme(:)
        logical, intent(in), optional :: start_with_B

        logical :: start_with_B_1
        logical :: opA      ! .true. for A, .false for B
        integer :: split_steps, k

        opA = .true.        
        if (present(start_with_B).and. start_with_B) then
            opA = .false.
        end if
       
        
        split_steps = ubound(splitting_scheme,1)
        do k=1, split_steps      
            if (opA) then
                call this%imaginary_time_propagate_A(dt*splitting_scheme(k))
            else
                if (k==1) then 
                    call this%imaginary_time_propagate_B(dt*splitting_scheme(k), 0)
                else if (k==split_steps) then 
                    call this%imaginary_time_propagate_B(dt*splitting_scheme(k), 1)
                else 
                    call this%imaginary_time_propagate_B(dt*splitting_scheme(k), 2)
                end if    
            end if
            opA = .not.opA
        end do
    end subroutine splitting_imaginary_time_step


    subroutine print_local_imaginary_time_orders(this, dt0, rows, &
                           start_with_B, splitting_scheme, extrapolation_order, fraunz)
       class(S(wf_schroedinger)), intent(inout) :: this
       real(kind=prec), intent(in) :: dt0
       integer, intent(in), optional :: rows 
       logical, intent(in), optional :: start_with_B
       _COMPLEX_OR_REAL_(kind=prec), intent(in), target, optional :: splitting_scheme(:)
       integer, intent(in), optional :: extrapolation_order 
       integer, intent(in), optional :: fraunz
  

       type(S(wf_schroedinger)) :: psi0, psi1, psi2, psi3
       integer :: row, k, fraunz1, rows1, extrapolation_order_1
       real(kind=prec) :: err, err_old, p
       _COMPLEX_OR_REAL_(kind=prec) :: dt1
       _COMPLEX_OR_REAL_(kind=prec), pointer :: scheme(:)
#ifdef _REAL_
       real(kind=prec), target :: default_scheme(3) =(/ 0.5_prec, 1.0_prec, 0.5_prec /)
#else
       complex(kind=prec), target :: default_scheme(3) = & 
                          (/ (0.5_prec,0.0_prec), (1.0_prec,0.0_prec), (0.5_prec,0.0_prec) /)
#endif        

       fraunz1 = 10
       if (present(fraunz)) then
           fraunz1 = fraunz
       end if    

       rows1 = 10
       if (present(rows)) then
           rows1 = rows
       end if         

       extrapolation_order_1 = 1
       if (present(extrapolation_order)) then
           extrapolation_order_1 = extrapolation_order
       end if    

       if (present(splitting_scheme)) then
            scheme => splitting_scheme
            extrapolation_order_1 = 1
       else     
            scheme => default_scheme
       endif

       select type (m=>this%m); class is (S(schroedinger))
       psi0 =  S(wf_schroedinger)(m)
       psi1 =  S(wf_schroedinger)(m)
       if (extrapolation_order_1>=2) then
          psi2 =  S(wf_schroedinger)(m)
          psi3 =  S(wf_schroedinger)(m)
       end if   
       end select
       
    
       call psi0%copy(this)
       call psi1%copy(psi0)
       dt1 = dt0
       write (*,*) "            dt         err      p"
       write (*,*) "----------------------------------"
       do row=1,rows1
           if (extrapolation_order_1==1) then
               do k=1,fraunz1
                   call psi1%splitting_imaginary_time_step(dt1/real(fraunz1, kind=prec), scheme, start_with_B)
               end do
               call this%splitting_imaginary_time_step(dt1, scheme, start_with_B)
           else
               do k=1,fraunz1
                   call psi1%extrapolation_imaginary_time_step(dt1/real(fraunz1, kind=prec), &
                                                             extrapolation_order_1, psi2, psi3, start_with_B)
               end do
               call this%extrapolation_imaginary_time_step(dt1, extrapolation_order_1, psi2, psi3, start_with_B)
           end if

           err = this%distance(psi1)
           if (row==1) then
               if (this_proc==0) then
                   write (*, '(I3,2E12.3)') row, real(dt1, prec), err
               end if    
           else
               p = log(err_old/err)/log(2.0_prec);
               if (this_proc==0) then
                   write (*, '(I3,2E12.3, F7.2)') row, real(dt1, prec), err, p
               end if    
           end if
           err_old = err;
           dt1 = dt1/2.0_prec
           call this%copy(psi0)
           call psi1%copy(psi0)
       end do

       !TODO finalize psi0, psi1, psi2, psi3

    end subroutine print_local_imaginary_time_orders




    subroutine compute_groundstate(this, dt0, tol, max_iters, &
                                                      start_with_B, splitting_scheme, extrapolation_order)
#ifdef _QUADPRECISION_
       use tssmq_common
#else
       use tssm_common
#endif       
       class(S(wf_schroedinger)), intent(inout) :: this
       real(kind=prec), intent(in) :: dt0 ! TODO should be optional
       real(kind=prec), intent(in) :: tol ! TODO should be optional
       integer, intent(in) :: max_iters   ! TODO should be optional
       logical, intent(in), optional :: start_with_B
       _COMPLEX_OR_REAL_(kind=prec), intent(in), target, optional :: splitting_scheme(:)
       integer, intent(in), optional :: extrapolation_order 

       type(S(wf_schroedinger)) :: psi0, psi1, psi2, psi3
       type(S(wf_schroedinger)) :: psi_old

       _COMPLEX_OR_REAL_(kind=prec), pointer :: scheme(:)
#ifdef _REAL_
       real(kind=prec), target :: default_scheme(3) =(/ 0.5_prec, 1.0_prec, 0.5_prec /)
#else
       complex(kind=prec), target :: default_scheme(3) = & 
                          (/ (0.5_prec,0.0_prec), (1.0_prec,0.0_prec), (0.5_prec,0.0_prec) /)
#endif                          
       integer :: k, k_check
       logical :: start_with_B_1
       integer :: extrapolation_order_1 
       _COMPLEX_OR_REAL_(kind=prec) :: dt
       real(kind=prec) :: E, E_mu, E_dev, err, E1, E_mu1, E_dev1, err1, N, E_est, E_old 
       real(kind=prec) :: ddd 
       integer(kind=8) :: calc
       real :: calc_time 

       call tick(calc)

       start_with_B_1 = .false.
       if (present(start_with_B)) then
           start_with_B_1 = start_with_B
       end if
       
       extrapolation_order_1 = 1
       if (present(extrapolation_order)) then
           extrapolation_order_1 = extrapolation_order
       end if    

       if (present(splitting_scheme)) then
            scheme => splitting_scheme
            extrapolation_order_1 = 1
       else     
            scheme => default_scheme
       endif


       select type (m=>this%m); class is (S(schroedinger))
       psi1 =  S(wf_schroedinger)(m)
       if (extrapolation_order_1>=2) then
          psi2 =  S(wf_schroedinger)(m)
          psi3 =  S(wf_schroedinger)(m)
       end if   
       psi_old =  S(wf_schroedinger)(m)       
       end select

#ifdef _REAL_
       dt = dt0
#else       
       dt = cmplx(dt0, 0.0_prec, prec)
#endif       
       ! initial value, better choices !??
       this%u = 1.0_prec
       this%is_real_space = .true.

       call psi_old%copy(this)

       k = 0
       k_check =20 
       E_old = 1e6_prec ! TODO infinity
       do k=0, max_iters
           if (mod(k,k_check)==0) then
               if (start_with_B_1) then
                   call this%to_real_space()
               else
                   call this%to_frequency_space()
               end if
               call psi1%copy(this)

               if (extrapolation_order_1==1) then
                   call psi1%splitting_imaginary_time_step(0.5_prec*dt, scheme, start_with_B_1)
               else
                   call psi1%extrapolation_imaginary_time_step(0.5_prec*dt, extrapolation_order_1, psi2, psi3, start_with_B_1)
                   !call psi1%extrapolation_imaginary_time_step(0.5_prec*dt, extrapolation_order_1, psi2, this, start_with_B_1)
               endif
#ifndef _REAL_
               call psi1%to_real_space()
               psi1%u = real(psi1%u, prec)
#endif 
               call psi1%normalize
               call psi1%get_energy_expectation_deviation(E_mu1, E_dev1)
               E1 = E_mu1 - psi1%interaction_energy()
               err1 = E_dev1/E1
           end if

           if (extrapolation_order_1==1) then
               call this%splitting_imaginary_time_step(dt, scheme, start_with_B_1)
           else
               call this%extrapolation_imaginary_time_step(dt, extrapolation_order_1, psi2, psi3, start_with_B_1)
               !call this%extrapolation_imaginary_time_step(dt, extrapolation_order_1, psi2, psi1, start_with_B_1)
           end if

#ifndef _REAL_
           call this%to_real_space()
           this%u = real(this%u, prec)
#endif 
!          call this%normalize(N)           
           call this%normalize           
           call this%get_energy_expectation_deviation(E_mu, E_dev)
           E = E_mu - this%interaction_energy()
           err = E_dev/E
!#ifdef _REAL_
!           E_est = -log(N)/dt
!#else
!           E_est = -log(N)/real(dt, prec)
!#endif           

           ddd = this%distance(psi_old)
           calc_time = tock(calc)
           if (mod(k,k_check)==0) then
               if (this_proc==0) then
                   !write (*, '(I5,2E24.15,2E12.3,E24.15,2E12.3)') k, E, E_mu,  E_old-E, err, E1, E_old-E1, err1
                   write (*, '(I5,2E24.15,3E12.3,F10.2,E24.15,2E12.3)') k, E, E_mu,  E_old-E, err, ddd, calc_time, &
                                                                       E1, E_old-E1, err1
               endif     
               !if (err1<err .or. E_old-E<1e-14_prec) then
               !if (err1<err .and. E1<E) then
               if (err1<err) then
               !if (E1<E) then
                   if (this_proc==0) then
#ifdef _REAL_
                       print *, "changed step size, old:", dt, "  new:", 0.5_prec*dt
#else
                       print *, "changed step size, old:", real(dt, prec), "  new:", 0.5_prec*real(dt, prec)
#endif
                   endif
                 
                   dt = 0.5_prec*dt
                   !if (extrapolation_order==1) then
                       call this%copy(psi1)
                       !Note that if extrapolation_step>1, then psi1 has been overwritten and
                       !cannot be used anymore
                   !end if    
               end if
           else    
               if (this_proc==0) then
                   !write (*, '(I5,2E24.15,2E12.3)') k, E, E_mu, E_old-E,  err
                   write (*, '(I5,2E24.15,3E12.3,F10.2)') k, E, E_mu, E_old-E,  err, ddd, calc_time 
               end if    
           end if    

           if (err<tol) exit

           E_old = E
           call psi_old%copy(this)
        
       end do

       !TODO finalize psi1, psi2, psi3

    end subroutine compute_groundstate

#ifdef _QUADPRECISION_
end module S(tssmq_schroedinger)
#else
end module S(tssm_schroedinger)
#endif
