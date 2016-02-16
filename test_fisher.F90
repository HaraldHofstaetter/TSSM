module tssm_fisher
    use tssm_fourier
    implicit none

    !------------------------------------------------
    !complex versions
    !------------------------------------------------

    type, extends(fourier_1D) :: fisher
        integer :: method_for_B = 0  
    end type fisher

    interface fisher
        module procedure :: new_fisher
    end interface fisher

    type, extends(wf_fourier_1D) :: wf_fisher
    contains 
        procedure :: clone => clone_wf_fisher
        procedure :: propagate_B => propagate_B_wf_fisher
    end type wf_fisher

    interface wf_fisher
        module procedure new_wf_fisher
    end interface wf_fisher

 
    !------------------------------------------------
    !real versions
    !------------------------------------------------
    type, extends(fourier_real_1D) :: fisher_real
    end type fisher_real

    interface fisher_real
        module procedure :: new_fisher_real
    end interface fisher_real

    type, extends(wf_fourier_real_1D) :: wf_fisher_real
    contains 
        procedure :: clone => clone_wf_fisher_real
        procedure :: propagate_B => propagate_B_wf_fisher_real
    end type wf_fisher_real

    interface wf_fisher_real
        module procedure new_wf_fisher_real
    end interface wf_fisher_real



contains 
    !------------------------------------------------
    !complex versions
    !------------------------------------------------

    function new_fisher(nx, xmin, xmax, boundary_conditions) result(this)
        type(fisher) :: this
        real(kind=prec), intent(in) :: xmin 
        real(kind=prec), intent(in) :: xmax
        integer, intent(in) :: nx
        integer, intent(in), optional :: boundary_conditions

        this%fourier_1D = fourier_1D(nx, xmin, xmax, boundary_conditions)
    end function new_fisher

    function new_wf_fisher(m) result(this)
        type(wf_fisher) :: this
        class(fisher), target, intent(inout) :: m
        
        this%wf_fourier_1D = wf_fourier_1D(m, coefficient=(1.0_prec, 0.0_prec))
    end function new_wf_fisher

    function clone_wf_fisher(this) result(clone)
        class(wf_fisher), intent(inout) :: this
        class(wave_function), pointer :: clone
        type(wf_fisher), pointer :: p

        allocate( p )
        select type (m=>this%m); class is (fisher)
        p = new_wf_fisher(m)
        end select
        clone => p
    end function clone_wf_fisher

    
    elemental function F(u)
        complex(kind=prec), intent(in) :: u
        complex(kind=prec) :: F
        F = u*(1.0_prec-u)
    end function F

    elemental function Phi(dt, u) result(u1)
        complex(kind=prec), intent(in) :: dt
        complex(kind=prec), intent(in) :: u
        complex(kind=prec) :: u1
        u1 = u/(u + exp(-dt)*(1.0_prec-u))
    end function Phi

    elemental function Heun_3rd_order(dt, u) result(u1)
        complex(kind=prec), intent(in) :: dt
        complex(kind=prec), intent(in) :: u
        complex(kind=prec) :: u1
        complex(kind=prec) :: k1, k2, k3 
        k1 = F(u)
        k2 = F(u+dt/3.0_prec*k1)
        k3 = F(u+dt*2.0_prec/3.0_prec*k2)
        u1 = u+dt*(0.25_prec*k1+0.75_prec*k3)
    end function Heun_3rd_order

    elemental function Runge_Kutta_4th_order(dt, u) result(u1)
        complex(kind=prec), intent(in) :: dt
        complex(kind=prec), intent(in) :: u
        complex(kind=prec) :: u1
        complex(kind=prec) :: k1, k2, k3, k4 
        k1 = F(u)
        k2 = F(u+dt*0.5_prec*k1)
        k3 = F(u+dt*0.5_prec*k2)
        k4 = F(u+dt*k3)
        u1 = u+dt*(1.0_prec/6.0_prec*k1+2.0_prec/6.0_prec*k2+2.0_prec/6.0_prec*k3+1.0_prec/6.0_prec*k4)
    end function Runge_Kutta_4th_order

    subroutine propagate_B_wf_fisher(this, dt)
        class(wf_fisher), intent(inout) :: this
        complex(kind=prec), intent(in) :: dt

       complex(kind=prec), dimension(:), pointer :: u

       select type (m=>this%m); class is (fisher)

       call this%to_real_space  

       u => this%u !(m%g%n1min:m%g%n1max)

       select case(m%method_for_B)
       case(1) ! explicit Euler
          u = u + dt*F(u)
       case(2) ! explicit midpoint rule  
          u = u + dt*F(u + 0.5_prec*dt*F(u)) 
       case(3)                                         
          u =  Heun_3rd_order(dt, u)
       case(4)
          u = Runge_Kutta_4th_order(dt, u)
       case default
          u = Phi(dt, u) ! exact
       end select   

       class default
           stop "E: wrong spectral method for fisher wave function"
       end select
    end subroutine propagate_B_wf_fisher

    function initial_value(x)
        complex(kind=prec) :: initial_value
        complex(kind=prec), intent(in) :: x

        initial_value = sin(2*pi*x)*exp(x) !*(1.0_prec, 1.0_prec)
    end function initial_value


    !------------------------------------------------
    !real versions
    !------------------------------------------------

    function new_fisher_real(nx, xmin, xmax, boundary_conditions) result(this)
        type(fisher_real) :: this
        real(kind=prec), intent(in) :: xmin 
        real(kind=prec), intent(in) :: xmax
        integer, intent(in) :: nx
        integer, intent(in), optional :: boundary_conditions

        this%fourier_real_1D = fourier_real_1D(nx, xmin, xmax, boundary_conditions)
    end function new_fisher_real

    function new_wf_fisher_real(m) result(this)
        type(wf_fisher_real) :: this
        class(fisher_real), target, intent(inout) :: m
        
        this%wf_fourier_real_1D = wf_fourier_real_1D(m, coefficient=1.0_prec)
    end function new_wf_fisher_real

    function clone_wf_fisher_real(this) result(clone)
        class(wf_fisher_real), intent(inout) :: this
        class(real_wave_function), pointer :: clone
        type(wf_fisher_real), pointer :: p

        allocate( p )
        select type (m=>this%m); class is (fisher_real)
        p = new_wf_fisher_real(m)
        end select
        clone => p
    end function clone_wf_fisher_real


    subroutine propagate_B_wf_fisher_real(this, dt)
        class(wf_fisher_real), intent(inout) :: this
        real(kind=prec), intent(in) :: dt
        integer :: jx
        real(kind=prec) :: e_dt

       real(kind=prec), dimension(:), pointer :: u

       select type (m=>this%m); class is (fisher_real)

       call this%to_real_space  

       u => this%u 
       u =u/(u + exp(-dt)*(1.0_prec-u))

       class default
           stop "E: wrong spectral method for fisher_real wave function"
       end select
    end subroutine propagate_B_wf_fisher_real


    function initial_value_real(x)
        real(kind=prec) :: initial_value_real
        real(kind=prec), intent(in) :: x

        initial_value_real = cos(pi*x)!*(1.0_prec, 1.0_prec)
    end function initial_value_real

end module tssm_fisher

program test_fisher
    !use, intrinsic :: IEEE_ARITHMETIC
    use tssm_fisher
    implicit none 

    type(fisher) :: method 
    type(wf_fisher) :: psi
    type(wf_fisher) :: psi_ex

    type(fisher_real) :: method_real 
    type(wf_fisher_real) :: psi_real
    type(wf_fisher_real) :: psi_ex_real

    real(prec) :: tend

    integer :: rows
    integer :: steps 
    complex(prec) :: scheme_A44c(7)
    complex(prec) :: scheme_A33c(5)
    complex(prec) :: scheme(6) = (/ &
               (0.201639688260407656_prec, 0.105972321241365172_prec),  &
               (0.387747410753696807_prec, 0.100071120693574555_prec),   &
               (0.410612900985895537_prec, -0.206043441934939727_prec), &
               (0.410612900985895537_prec, -0.206043441934939727_prec), &
               (0.387747410753696807_prec, 0.100071120693574555_prec),  &
               (0.201639688260407656_prec, 0.105972321241365172_prec)   /)

    complex(prec) :: omega1, omega2

    omega1 = 1.0_prec/(2.0_prec - 2.0_prec**(1.0_prec/3.0_prec)*exp(cmplx(0.0_prec, 2.0_prec*pi/3.0_prec, prec)))
    omega2 = 1.0_prec -2.0_prec*omega1
    scheme_A44c = (/ 0.5_prec*omega1, omega1, 0.5_prec*(omega1+omega2), omega2, 0.5_prec*(omega1+omega2), omega1, 0.5_prec*omega1 /)  
    scheme_A33c = (/ cmplx(0.25_prec, sqrt(3.0_prec)/12.0_prec, prec), &
                     cmplx(0.5_prec, sqrt(3.0_prec)/6.0_prec, prec),   &
                     cmplx(0.5_prec, 0.0_prec, prec), &
                     cmplx(0.5_prec, -sqrt(3.0_prec)/6.0_prec, prec),  &
                     cmplx(0.25_prec, -sqrt(3.0_prec)/12.0_prec, prec) /)


    !call ieee_set_underflow_mode(.true.)
#if 0
!----------------------
print *, "REAL"
    method_real = fisher_real(nx = 128,          &
                    xmin = 0.0_prec,  &
                    xmax = 1.0_prec,  &
                    boundary_conditions = dirichlet)
    psi_real =  wf_fisher_real(method_real)
    psi_ex_real =  wf_fisher_real(method_real)
    
    tend = 0.1_prec

    rows =  10
    steps = 10*2**rows
    call psi_ex_real%set(initial_value_real)
    call psi_ex_real%run(dt = tend/steps,& 
                    t0 = 0.0_prec, &
                    steps=steps, &
                    splitting_scheme = (/ 0.5_prec, 1.0_prec, 0.5_prec /), &
                    operator_sequence="AB")


    call psi_real%set(initial_value_real)
    call psi_real%print_orders(reference_solution=psi_ex_real, tend=tend, rows=rows, dt=tend, &
            splitting_scheme = (/ 0.5_prec, 1.0_prec, 0.5_prec /), &
      !     complex_splitting_scheme = (/ (8.0_prec,-1.0_prec)/26.0_prec, (18.0_prec,-1.0_prec)/25, &
       !                                   (18.0_prec, 1.0_prec)/26, (7.0_prec, 1.0_prec)/25 /), &
                    operator_sequence="AB")

#endif
!----------------------
print *
print *, "COMPLEX"

    method = fisher(nx = 1024,          &
                    xmin = 0.0_prec,  &
                    xmax = 1.0_prec,  &
                    boundary_conditions = dirichlet)

    psi =  wf_fisher(method)
    psi_ex =  wf_fisher(method)

#if 1
    rows = 20
    tend = 0.01_prec
    call psi%set(initial_value)
    method%method_for_B = 0
    call psi%print_local_orders(rows=rows, dt=tend, &
!            splitting_scheme = (/ 0.5_prec, 1.0_prec, 0.5_prec /), &
            !splitting_scheme = (/ 1.0_prec, 1.0_prec /), &
!           complex_splitting_scheme = (/ (8.0_prec,-1.0_prec)/26.0_prec, (18.0_prec,-1.0_prec)/25, &
!                                          (18.0_prec, 1.0_prec)/26, (7.0_prec, 1.0_prec)/25 /), &
                    complex_splitting_scheme = scheme_A44c, &
                    operator_sequence="AB")

#endif


#if 0
    tend = 0.1_prec

    rows = 10
    steps = 10*2**rows
    call psi_ex%set(initial_value)
    method%method_for_B = 0
    call psi_ex%run(dt = tend/steps,&
                    t0 = 0.0_prec, &
                    steps=steps, &
!                    splitting_scheme = (/ 0.5_prec, 1.0_prec, 0.5_prec /), &
!           complex_splitting_scheme = (/ (8.0_prec,-1.0_prec)/26.0_prec, (18.0_prec,-1.0_prec)/25, &
!                                          (18.0_prec, 1.0_prec)/26, (7.0_prec, 1.0_prec)/25 /), &
                    complex_splitting_scheme = scheme_A44c, &
                    operator_sequence="AB")


    tend = 0.1_prec
    rows = 15
    call psi%set(initial_value)
    method%method_for_B = 0
    call psi%print_orders(reference_solution=psi_ex, tend=tend, rows=rows, dt=tend, &
!            splitting_scheme = (/ 0.5_prec, 1.0_prec, 0.5_prec /), &
            !splitting_scheme = (/ 1.0_prec, 1.0_prec /), &
!            splitting_scheme = (/ 0.5_prec, 1.0_prec, 0.5_prec /), &
!           complex_splitting_scheme = (/ (8.0_prec,-1.0_prec)/26.0_prec, (18.0_prec,-1.0_prec)/25, &
!                                          (18.0_prec, 1.0_prec)/26, (7.0_prec, 1.0_prec)/25 /), &
                    complex_splitting_scheme = scheme_A44c, &
                    operator_sequence="AB")
#endif


end program test_fisher
