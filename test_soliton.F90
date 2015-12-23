module soliton_module
    use tssm, only: prec
    implicit none 

    real(prec) :: a = 2_prec
    real(prec) :: b = 1_prec;
    real(prec) :: c = 0_prec;

contains
    function soliton(x, t)
        complex(kind=prec) :: soliton
        real(kind=prec), intent(in) :: x
        real(kind=prec), intent(in) :: t
        real(kind=prec) :: h
        h = .5_prec*(a**2 - b**2)*t - b*x
        soliton = (a/cosh(a*(b*t+x-c))) * exp(cmplx(0_prec, h, prec)) 
    end function soliton
end module soliton_module



program test_soliton
    use tssm_schroedinger
    use soliton_module  
    use tssm_splitting_schemes

    implicit none 

    type(schroedinger_1D) :: method 
    type(wf_schroedinger_1D) :: psi
    type(wf_schroedinger_1D) :: psi_ex
    real(prec) :: tend, dt 

    call initialize_tssm

    method = schroedinger_1D(nx = 1024,          &
                       xmin = -16.0_prec,  &
                       xmax =  16.0_prec,  &
                       hbar = 1.0_prec,    &
                       mass = 1.0_prec,    &
                       cubic_coupling = -1.0_prec, &
                       boundary_conditions = dirichlet)
    psi =  wf_schroedinger_1D(method)
    psi_ex =  wf_schroedinger_1D(method)

    
    tend = +1.0_prec

    !get initial solution:
    call psi%set_t(soliton, 0.0_prec) 

    !get exact final solution
    call psi_ex%set_t(soliton, tend) 
   
    dt =  1.0e-1_prec
    call psi%print_orders(reference_solution=psi_ex, t0=0.0_prec, tend=tend, rows=12, dt=tend, &
!            splitting_scheme = (/ 0.5_prec, 1.0_prec, 0.5_prec /), &
            splitting_scheme = coeffs_15A, &
            operator_sequence="AB")
!    call psi%save("xxx.h5")

!    dt =  1.0e-1_prec
!    call psi%run_adaptive(dt, (/ 0.0_prec, tend /), &
!                    tol = 1e-9_prec, &
!                    !splitting_scheme = coeffs_ABC_6, & !coeffs_17, & ! palindromic
!                    splitting_scheme = coeffs_15, &
!                    associated_scheme = coeffs_15A, &
!                    order = 4, &
!                    operator_sequence = "AB", &
!                    solution_out = solution_out_schroedinger_1D)
!    print *, "ERROR", psi%distance(psi_ex)


    call finalize_tssm
                  
end program test_soliton
