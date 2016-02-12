module harmonic_module
    use tssm_schroedinger
    implicit none 

    real(prec), parameter :: omega_x = 1.1_prec
    real(prec), parameter :: hbar = 1.1_prec
    real(prec), parameter :: mass = 1.1_prec
    real(prec), parameter :: x0 = sqrt(hbar/mass/omega_x)
    real(prec), parameter :: alpha = 1.1_prec

    real(prec) :: t = 0.0_prec

contains
    function harmonic_potential(x)
        real(kind=prec) :: harmonic_potential
        real(kind=prec), intent(in) :: x
        harmonic_potential = .5_prec*mass*omega_x**2*x**2;
    end function harmonic_potential

    function coherent_state(x)
        complex(kind=prec) :: coherent_state
        real(kind=prec), intent(in) :: x
        
        real(kind=prec) :: xi 
        complex(kind=prec) :: alpha_t

        xi = x/x0
        alpha_t=alpha*exp(cmplx(0.0_prec, -omega_x*t, kind=prec))

        coherent_state =  1.0_prec/sqrt(x0*sqrt(pi))*exp(cmplx(0.0_prec, -0.5_prec*omega_x*t, kind=prec)) &
                          *exp(sqrt(2.0_prec)*alpha_t*xi - 0.5_prec*xi**2 - real(alpha_t, kind=prec)*alpha_t)
    end function coherent_state

    function observable_x(x)
        real(kind=prec) :: observable_x
        real(kind=prec), intent(in) :: x
        observable_x = x
    end function 
end module harmonic_module



program test_harmonic
    use harmonic_module
    implicit none 

    type(schroedinger_1D) :: method 
    type(wf_schroedinger_1D) :: psi
    type(wf_schroedinger_1D) :: psi_ex
    real(prec) :: tend, E_kin, E_pot

    call initialize_tssm

    method = schroedinger_1D(nx = 256,          &
                       xmin = -8.0_prec,  &
                       xmax =  8.0_prec,  &
                       hbar = hbar,   &
                       mass = mass,   &
                       potential = harmonic_potential, &
                       boundary_conditions = periodic)

    psi =  wf_schroedinger_1D(method)
    psi_ex =  wf_schroedinger_1D(method)

    
    tend = 1.0_prec

    !get initial solution:
    t = 0_prec
    call psi%set(coherent_state) 
    E_kin = psi%kinetic_energy()
    E_pot = psi%potential_energy()
    print *, "E_kin", E_kin 
    print *, "E_pot", E_pot
    print *, "E_tot", E_kin+E_pot
    print *, "norm", psi%norm()
    !call psi%save("xxx0.h5")

    !get exact final solution
    t = tend
    call psi_ex%set(coherent_state) 
    E_kin = psi_ex%kinetic_energy()
    E_pot = psi_ex%potential_energy()
    print *, "E_kin", E_kin 
    print *, "E_pot", E_pot
    print *, "E_tot", E_kin+E_pot
    print *, "norm", psi_ex%norm()
    print *, "observable V", psi_ex%observable(harmonic_potential)
    print *, "observable x", psi_ex%observable(observable_x)
    !call psi_ex%save("xxx_ex.h5")

 !   call psi%run(dt=tend/1000.0_prec, t0= 0.0_prec, splitting_scheme = (/ 0.5_prec, 1.0_prec, 0.5_prec /), &
 !                 steps=1000, start_with_B=.false., solution_out = solution_out_schroedinger_1D)
   ! call psi%save("xxx.h5")
!    print *, "err", psi%distance(psi_ex)

 
    call psi%print_orders(reference_solution=psi_ex, t0=0.0_prec, tend=tend, rows=10, dt=tend, &
           splitting_scheme = (/ 1.0_prec, -1.0_prec/24, -2.0_prec/3, 3.0_prec/4, 2.0_prec/3, 7.0_prec/24 /), &
!            complex_splitting_scheme = (/ (8.0_prec,-1.0_prec)/26.0_prec, (18.0_prec,-1.0_prec)/25, &
!                                          (18.0_prec, 1.0_prec)/26, (7.0_prec, 1.0_prec)/25 /), &
          operator_sequence="AB")

                  
    call finalize_tssm
   
end program test_harmonic
