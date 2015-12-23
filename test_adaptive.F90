program test_adaptive
    use tssm_schroedinger
    implicit none

    real(prec), parameter :: cubic_coupling = 390.0_prec
    !real(prec), parameter :: cubic_coupling = 406.89_prec

    call initialize_tssm

    call compute_groundstate
    call run
    
    call finalize_tssm

contains
    subroutine compute_groundstate
        type(schroedinger_real_1D) :: method 
        type(wf_schroedinger_real_1D) :: psi

        method = schroedinger_real_1D(nx = 2048, &
                       xmin = -100.0_prec,       &
                       xmax =  100.0_prec,       &
                       hbar = 1.0_prec,          &
                       mass = 1.0_prec,          &
                       cubic_coupling = cubic_coupling, &
                       potential = harmonic_trap_1D,    &
                       boundary_conditions = periodic)

       psi =  wf_schroedinger_real_1D(method)     

       call psi%compute_groundstate(dt0 = 0.05_prec, &
                       tol = 1e-10_prec,    &
                       max_iters = 10000,   &
                       extrapolation_order = 2)

       call psi%save("groundstate.h5")                       

       call psi%finalize
       
    end subroutine compute_groundstate

    
    subroutine run
        use tssm_splitting_schemes
        type(schroedinger_1D) :: method 
        type(wf_schroedinger_1D) :: psi, psi0

        real(prec) :: dt, t0, tend        


        method = schroedinger_1D(nx = 4*2048, &
                       xmin = -400.0_prec,       &
                       xmax =  400.0_prec,       &
                       hbar = 1.0_prec,          &
                       mass = 1.0_prec,          &
                       cubic_coupling = cubic_coupling, &
                       potential = weak_periodic_1D,    &
                       boundary_conditions = periodic)

       psi =  wf_schroedinger_1D(method)     
       psi0 =  wf_schroedinger_1D(method)     

       call psi%load("groundstate.h5")  
       call psi0%copy(psi)

!       call psi%print_local_orders(rows=10, dt=0.1_prec, splitting_scheme = coeffs_15A)

       t0 = 0.0_prec
       tend = 1_prec
       dt = 1.0e-2_prec


       call psi%run_adaptive(dt, (/ t0, tend /), &
                    tol = 1e-7_prec, &
                    !splitting_scheme = coeffs_17, & ! palindromic
                    splitting_scheme = coeffs_15A, &
                    controller_scheme = coeffs_15, &
                    order = 4, &
                    operator_sequence = "BA", &
                    solution_out = solution_out_schroedinger_1D)

!       call psi%save("solution_t30.h5")       
!       call psi%load("solution_t30.h5")       

       call psi%run_adaptive(dt, (/ tend, t0 /), &
                    tol = 1e-7_prec, &
                    !splitting_scheme = coeffs_17, & ! palindromic
                    splitting_scheme = coeffs_15A, &
                    controller_scheme = coeffs_15, &
                    order = 4, &
                    operator_sequence = "BA", &
                    solution_out = solution_out_schroedinger_1D)

 !      call psi%save("solution_t0_back.h5")       
print *, "ERROR", psi%distance(psi0)



       call psi%finalize

        
    end subroutine run



    function harmonic_trap_1D(x) result(V)
        real(kind=prec), intent(in) :: x
        real(kind=prec) :: V
        V = 0.5_prec * x**2
    end function harmonic_trap_1D

    function weak_periodic_1D(x) result(y)
        real(kind=prec), intent(in) :: x
        real(kind=prec) :: y
        y = 1.4_prec * cos(11.46_prec*x) 
    end function weak_periodic_1D

    

end program test_adaptive
