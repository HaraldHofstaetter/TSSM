program test_hermite
    !use, intrinsic :: IEEE_ARITHMETIC
    use tssm_hermite
    use tssm_fourier
    use tssm_schroedinger
    implicit none 
    
    real(prec) :: t

#if 0

    type(hermite_1D) :: method 
    type(fourier_1D) :: method_f

    type(wf_hermite_1D) :: psi, psi_ex
    type(wf_fourier_1D) :: psi_f
 
    method = hermite_1D(100, 1.0_prec)
    method_f = fourier_1D(100, xmin=-20.0_prec,xmax=+20.0_prec)

    psi = wf_hermite_1D(method, coefficient=(0_prec, 1.0_prec))
    psi_ex = wf_hermite_1D(method)
    psi_f = wf_fourier_1D(method_f)

    t = 0.0_prec

    call psi%set(coherent_state)
    call psi_ex%set(coherent_state)
    call psi_f%set(coherent_state)

    call psi%save("xxx.h5")
    call psi_f%save("xxx_f.h5")

    t = 1.1_prec
    call psi%propagate_A(cmplx(t,kind=prec))
    call psi%save("xxx1.h5")
    call psi_ex%set(coherent_state)
    call psi_ex%save("xxx2.h5")
    print *, "N", psi%norm(), psi_ex%norm()

    print *, "ERR", psi%distance(psi_ex)
#else
    real(prec) :: mu, E_dev 
    type(schroedinger_hermite_3D) :: method 
    type(schroedinger_3D) :: method_f

    type(wf_schroedinger_hermite_3D) :: psi
    type(wf_schroedinger_3D) :: psi_f
 
    method = schroedinger_hermite_3D(52, 1.1_prec, &
                                     54,  1.2_prec, &
                                     56,  1.3_prec, &
                                     hbar =1.4_prec, &
                                     mass =1.5_prec, &
                                     cubic_coupling = 0001.0_prec)
    method_f = schroedinger_3D(nx=128, xmin=-20.0_prec, xmax=+20.0_prec, &
                               ny=128, ymin=-20.0_prec, ymax=+20.0_prec, &
                               nz=128, zmin=-20.0_prec, zmax=+20.0_prec, &
                               hbar=1.4_prec, &
                               mass=1.5_prec, &
                               potential = harmonic_3D, &
                               boundary_conditions = periodic, &
                               cubic_coupling = 0001.0_prec)

    psi = wf_schroedinger_hermite_3D(method)
    psi_f = wf_schroedinger_3D(method_f)


    call psi%rset(f_3D)
    call psi_f%rset(f_3D)


    print *, "E", psi%kinetic_energy()+ psi%potential_energy()+ psi%interaction_energy(), &
                  psi_f%kinetic_energy() + psi_f%potential_energy()+ psi_f%interaction_energy()

    call psi%get_energy_expectation_deviation(mu, E_dev)
    print *, "E1", mu-psi%interaction_energy()
    call psi_f%get_energy_expectation_deviation(mu, E_dev)
    print *, "E1", mu-psi_f%interaction_energy()
    print *, "E_kin", psi_f%kinetic_energy()
    print *, "E_pot", psi_f%potential_energy()
    print *, "E_int", psi_f%interaction_energy()
    print *, "norm", psi_f%norm()

    



#endif

contains
    function f_1D(x) result(y)
       real(kind=prec), intent(in) :: x
       real(kind=prec) :: y

       y = exp(-x**2)
    end function f_1D

    function harmonic_1D(x) result(y)
       real(kind=prec), intent(in) :: x
       real(kind=prec) :: y

       y = 0.5_prec*1.2_prec**2*x**2
    end function harmonic_1D

   function harmonic_2D(x,y) result(z)
       real(kind=prec), intent(in) :: x,y
       real(kind=prec) :: z

       z = 0.5_prec*(1.2_prec**2*x**2 + 1.4_prec**2*y**2)
    end function harmonic_2D

    function f_2D(x,y) result(z)
       real(kind=prec), intent(in) :: x,y
       real(kind=prec) :: z

       z = exp(-.7_prec*x**2-.9_prec*y**2)
    end function f_2D

   function harmonic_3D(x,y,z) result(w)
       real(kind=prec), intent(in) :: x,y,z
       real(kind=prec) :: w

       w = 0.5_prec*(1.1_prec**2*x**2 + 1.2_prec**2*y**2 + 1.3_prec**2*z**2)
    end function harmonic_3D

    function f_3D(x,y,z) result(w)
       real(kind=prec), intent(in) :: x,y,z
       real(kind=prec) :: w

       w = exp(-.8_prec*x**2-.7_prec*y**2-.9*z**2)
    end function f_3D



    function coherent_state(x)
        complex(kind=prec) :: coherent_state
        real(kind=prec), intent(in) :: x
        
        real(kind=prec) :: xi 
        complex(kind=prec) :: alpha_t 
        
        real(prec), parameter :: alpha = 1.3_prec
        real(prec), parameter :: omega_x = 1.0_prec
        real(prec), parameter :: x0 = sqrt(1.0_prec/omega_x)

        xi = x/x0
        alpha_t=alpha*exp(cmplx(0.0_prec, -omega_x*t, kind=prec))

        coherent_state =  1.0_prec/sqrt(x0*sqrt(pi))*exp(cmplx(0.0_prec, -0.5_prec*omega_x*t, kind=prec)) &
                          *exp(sqrt(2.0_prec)*alpha_t*xi - 0.5_prec*xi**2 - real(alpha_t, kind=prec)*alpha_t)
    end function coherent_state


end program test_hermite

