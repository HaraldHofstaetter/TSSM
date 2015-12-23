module tssm_coupled_nonlin_schroe
    use tssm_fourier
    use tssm_schroedinger
    implicit none

    type, extends(multicomponent_fourier_1D) :: coupled_nonlin_schroe
        real(kind=prec) :: g11
        real(kind=prec) :: g12
        real(kind=prec) :: g21
        real(kind=prec) :: g22
    end type coupled_nonlin_schroe

    interface coupled_nonlin_schroe
        module procedure :: new_coupled_nonlin_schroe
    end interface coupled_nonlin_schroe

    type, extends(wf_multicomponent_fourier_1D) :: wf_coupled_nonlin_schroe
    contains 
        procedure :: clone => clone_wf_coupled_nonlin_schroe
        procedure :: propagate_B => propagate_B_wf_coupled_nonlin_schroe
        !procedure :: propagate_C => propagate_C_wf_coupled_nonlin_schroe
    end type wf_coupled_nonlin_schroe

    interface wf_coupled_nonlin_schroe
        module procedure new_wf_coupled_nonlin_schroe
    end interface wf_coupled_nonlin_schroe
contains

    function new_coupled_nonlin_schroe(nx, xmin, xmax, &
                                      g11, g12, g21, g22) result(this)
        type(coupled_nonlin_schroe) :: this
        integer, intent(in) :: nx
        real(kind=prec), intent(in) :: xmin 
        real(kind=prec), intent(in) :: xmax
        real(kind=prec), intent(in) :: g11, g12, g21, g22 

        type(schroedinger_1D), pointer :: b
        this%g11 = g11
        this%g12 = g12
        this%g21 = g21
        this%g22 = g22
        this%multicomponent_fourier_1D = new_multicomponent_fourier_1D(  &
                           (/ (1.0_prec, 0.0_prec),  (1.0_prec, 0.0_prec) /), &
                           nx, xmin, xmax, periodic)
        allocate( b )
        b = schroedinger_1D(nx, xmin, xmax, &
               hbar = 1.0_prec, &
               mass = 1.0_prec, &
               boundary_conditions = periodic)
        if (.not.associated(b%V)) then
             call b%g%allocate_real_gridfun(b%V)
        end if    
        this%b => b
        this%dset_names_real(1) = "psi1_real"
        this%dset_names_imag(1) = "psi1_imag"
        this%dset_names_real(2) = "psi2_real"
        this%dset_names_imag(2) = "psi2_imag"

    end function new_coupled_nonlin_schroe

    function new_wf_coupled_nonlin_schroe(m) result(this)
        type(wf_coupled_nonlin_schroe) :: this
        class(coupled_nonlin_schroe), target, intent(inout) :: m
        integer :: j
        type(wf_schroedinger_1D), pointer :: p
        
        this%m => m
        allocate( this%C( m%nc ) )
        select type (b=>m%b)
        class is (schroedinger_1D)
        do j = 1, m%nc
            allocate(  p )
            p = new_wf_schroedinger_1D(b)
            this%C(j)%p => p
        end do    
        end select
    end function new_wf_coupled_nonlin_schroe

   function clone_wf_coupled_nonlin_schroe(this) result(clone)
        class(wf_coupled_nonlin_schroe), intent(inout) :: this
        class(wave_function), pointer :: clone
        type(wf_coupled_nonlin_schroe), pointer :: p
        
        select type (m=>this%m); class is (coupled_nonlin_schroe)
        allocate( p )
        p = wf_coupled_nonlin_schroe(m)
        end select
        clone => p
    end function clone_wf_coupled_nonlin_schroe


    subroutine propagate_B_wf_coupled_nonlin_schroe(this, dt)
        class(wf_coupled_nonlin_schroe), intent(inout) :: this
        complex(kind=prec), intent(in) :: dt

        complex(kind=prec), pointer :: u1(:)
        complex(kind=prec), pointer :: u2(:)
#ifdef _OPENMP
        integer :: j
#endif     
        select type (m=>this%m); class is (coupled_nonlin_schroe)
        select type (b=>m%b); class is (schroedinger_1D)
        call this%to_real_space()

#ifndef _OPENMP
        u1 => this%C(1)%p%u
        u2 => this%C(2)%p%u
        u1 = exp((-dt*(0.0_prec, 1.0_prec))*(m%g11*( real(u1, prec)**2+aimag(u1)**2 ) &
                                            +m%g12*( real(u2, prec)**2+aimag(u2)**2 ))) * u1
        u2 = exp((-dt*(0.0_prec, 1.0_prec))*(m%g21*( real(u1, prec)**2+aimag(u1)**2 ) &
                                            +m%g22*( real(u2, prec)**2+aimag(u2)**2 ))) * u2
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
        do j=1,n_threads
            u1 => this%C(1)%p%u(lbound(this%C(1)%p%u,1)+b%g%jj(j-1):lbound(this%C(1)%p%u,1)+b%g%jj(j)-1)
            u2 => this%C(2)%p%u(lbound(this%C(2)%p%u,1)+b%g%jj(j-1):lbound(this%C(2)%p%u,1)+b%g%jj(j)-1)
            u1 = exp((-dt*(0.0_prec, 1.0_prec))*(m%g11*( real(u1, prec)**2+aimag(u1)**2 ) &
                                                +m%g12*( real(u2, prec)**2+aimag(u2)**2 ))) * u1
            u2 = exp((-dt*(0.0_prec, 1.0_prec))*(m%g21*( real(u1, prec)**2+aimag(u1)**2 ) &
                                                +m%g22*( real(u2, prec)**2+aimag(u2)**2 ))) * u2
        end do    
!$OMP END PARALLEL DO
#endif

        end select
        end select
    end subroutine propagate_B_wf_coupled_nonlin_schroe




    
end module tssm_coupled_nonlin_schroe


program test_coupled_nonlin_schroe
    use tssm_coupled_nonlin_schroe
    use tssm_splitting_schemes
    implicit none

    ! parameters for exact solution (global varibles)
    real(kind=prec) :: t, delta, e, alpha, v
    real(kind=prec) :: x1, x2, v1, v2, alpha1, alpha2 

    call initialize_tssm

    call local_order_tables
    call order_tables
    !call crossing_solitons 
    
    call finalize_tssm

contains

subroutine order_tables

    real(prec) :: t0, tend, dt

    type(coupled_nonlin_schroe) :: method 
    type(wf_coupled_nonlin_schroe) :: psi, psi_ex     


   delta = 0.0_prec
   v = 1.1_prec
   e = 0.8_prec
   alpha = 1.0_prec
    
    method = coupled_nonlin_schroe(nx = 1024, &
                    xmin = -50.0_prec,  &
                    xmax = 50.0_prec,  &
                    g11 = -1.0_prec, g12 = -e, &
                    g21 = -e       , g22 = -1.0_prec) 

   psi =  wf_coupled_nonlin_schroe(method)
   psi_ex =  wf_coupled_nonlin_schroe(method)

   t0 = 0.0_prec
   tend = 5.0_prec
   dt = 0.1_prec

   t = t0
   call psi%set(exact_solution)
   call psi%save("xxx0.h5")

   t = tend
   call psi_ex%set(exact_solution)
   call psi_ex%save("xxx_ex.h5")

   call psi%print_orders(psi_ex, t0, tend, 8, dt, &
                         splitting_scheme = coeffs_13, &
                         operator_sequence = "AB") 

         
end subroutine order_tables


subroutine local_order_tables

    real(prec) :: t0, dt

    type(coupled_nonlin_schroe) :: method 
    type(wf_coupled_nonlin_schroe) :: psi


    delta = 0.0_prec
    v = 1.1_prec
    e = 0.8_prec
    alpha = 1.0_prec
    
    method = coupled_nonlin_schroe(nx = 1024, &
                    xmin = -50.0_prec,  &
                    xmax = 50.0_prec,  &
                    g11 = -1.0_prec, g12 = -e, &
                    g21 = -e       , g22 = -1.0_prec) 

   psi =  wf_coupled_nonlin_schroe(method)

   t0 = 0.0_prec
   dt = 0.1_prec
   call psi%set(exact_solution)

   call psi%print_local_orders(rows=8, dt=dt, &
                         splitting_scheme = coeffs_13, &
                         !controller_scheme = coeffs_9A, &
                         controller_scheme = palindromic, &
                         operator_sequence = "AB", &
                         get_reference_solution = get_exact_solution) 

end subroutine local_order_tables

subroutine get_exact_solution(psi, t1)
    class(wave_function), intent(inout) :: psi
    real(kind=prec), intent(in) :: t1 
    ! dummy variable called 't1' because we need access to global variable 't'    
    
    select type (psi); class is (wf_coupled_nonlin_schroe)    
    t = t1
    call psi%set(exact_solution)
    end select
end subroutine get_exact_solution


subroutine crossing_solitons 
    type(coupled_nonlin_schroe) :: method 
    type(wf_coupled_nonlin_schroe) :: psi

    real(prec) :: t0, tend, dt

    delta = 0.0_prec
    e = 0.7_prec ! 2.0_prec/3.0_prec

    method = coupled_nonlin_schroe(nx = 8192, &
                    xmin = -20.0_prec,  &
                    xmax = 80.0_prec,  &
                    g11 = -1.0_prec, g12 = -e, &
                    g21 = -e       , g22 = -1.0_prec) 

    psi =  wf_coupled_nonlin_schroe(method)

    t0 = 0.0_prec
    tend = 10.0_prec
    dt = 0.05_prec

    x1 = 0.0_prec
    v1 = 1.0_prec
    alpha1 = 0.5_prec
    x2 = +5.0_prec
    v2 = 0.1_prec
    alpha2 = 1.0_prec

    !t = t0
    !v = v1
    !alpha = alpha1

    call psi%set(initial_value)
    !call psi%set(exact_solution)
    call psi%save("xxx0.h5")

    call psi%run_adaptive(times = (/ t0, tend /), dt = dt, &
    !call psi%run(times = (/ t0, tend /), dt = dt, &
                tol = 1.0e-7_prec, &
                splitting_scheme = coeffs_15A, &
                controller_scheme = coeffs_15, &
                order = 4, &
                operator_sequence = "AB", &
                solution_out = solution_out_adaptive)
    call psi%save("xxx.h5")
end subroutine crossing_solitons 



   function exact_solution(x, which_component) result(psi)
        complex(kind=prec) :: psi
        real(kind=prec), intent(in) :: x
        integer, intent(in) :: which_component

        real(kind=prec) :: h1, h2

        h1 =sqrt(2.0_prec*alpha/(1.0_prec+e))/cosh(sqrt(2.0_prec*alpha)*(x-v*t))
        h2 = 0.5_prec*(v**2-delta**2)-alpha

        if (which_component==1) then
           psi = h1*exp(cmplx(0.0_prec, (v-delta)*x - h2*t, prec))
        else
           psi = h1*exp(cmplx(0.0_prec, (v+delta)*x - h2*t, prec))
        end if   
    end function exact_solution

   function initial_value(x, which_component) result(psi)
        complex(kind=prec) :: psi
        real(kind=prec), intent(in) :: x
        integer, intent(in) :: which_component
        real(kind=prec) :: h1, h2

        h1 =sqrt(2.0_prec*alpha1/(1.0_prec+e))/cosh(sqrt(2.0_prec*alpha1)*(x-x1))
        h2 =sqrt(2.0_prec*alpha2/(1.0_prec+e))/cosh(sqrt(2.0_prec*alpha2)*(x-x2))

        if (which_component==1) then
           psi = h1*exp(cmplx(0.0_prec, (v1-delta)*(x-x1), prec)) &
                +h2*exp(cmplx(0.0_prec, (v2-delta)*(x-x2), prec)) 
        else
           psi = h1*exp(cmplx(0.0_prec, (v1+delta)*(x-x1), prec)) &
                +h2*exp(cmplx(0.0_prec, (v2+delta)*(x-x2), prec)) 
        end if   
    end function initial_value 


    subroutine solution_out_adaptive(psi, t, dt, flag, abort_flag)
        class(wave_function), intent(inout) :: psi
        real(kind=prec), intent(in) :: t 
        real(kind=prec), intent(in) :: dt 
        character, intent(in) :: flag
        logical, intent(out) :: abort_flag
        integer, save :: step

        integer(kind=8), save :: calc 
        real :: calc_time         
        
        select type (psi); class is (wf_coupled_nonlin_schroe)

        if (flag=="i" .or. flag=="I") then ! init =======================
            call tick(calc)
            step = 0
            if (this_proc==0) then
                write (*,*), "    #                   t                  dt                cumm.time"
                write (*,*), "-----------------------------------------------------------------------" 
            end if           
        end if          
        
        if (flag=="i" .or. flag=="I" .or. flag=="e" .or. flag=="E") then ! execute =====
         calc_time = tock(calc)
            if (this_proc==0) then
              write (*, '(I6,2E20.12,F10.2)') step, t,  dt, calc_time 
          end if
            step = step + 1
        end if
    
        if (flag=="f" .or. flag=="F") then ! finish ====================
        end if   

        end select

        abort_flag = .false.
    end subroutine solution_out_adaptive


end program test_coupled_nonlin_schroe

