module gross_pitaevskii_module
    use tssm_schroedinger
    implicit none

contains
    function harmonic_trap_1D(x) result(V)
        real(kind=prec), intent(in) :: x
        real(kind=prec) :: V
        V = 0.5_prec * x**2
    end function harmonic_trap_1D

    function harmonic_trap_2D(x, y) result(V)
        real(kind=prec), intent(in) :: x, y
        real(kind=prec) :: V
        V = 0.5_prec * (x**2 + y**2)
    end function harmonic_trap_2D

    function weak_periodic(x) result(y)
        real(kind=prec), intent(in) :: x
        real(kind=prec) :: y
        y = 1.4_prec * cos(11.46_prec*x) 
    end function weak_periodic

    function init(x) result(y)
        real(kind=prec), intent(in) :: x
        real(kind=prec) :: y
        y =  exp(-x**2)
    end function init

end module gross_pitaevskii_module

program test_gross_pitaevskii
    use gross_pitaevskii_module
    implicit none

#if 1
    type(schroedinger_real_2D) :: method 
    type(wf_schroedinger_real_2D) :: psi

!    type(schroedinger_hermite_real_2D) :: method_h 
!    type(wf_schroedinger_hermite_real_2D) :: psi_h

    call initialize_tssm


    method = schroedinger_real_2D(nx = 256,     &
                       xmin = -20.0_prec,  &
                       xmax =  20.0_prec,  &
                       ny = 256, &
                       ymin = -20.0_prec,  &
                       ymax =  20.0_prec,  &
                       hbar = 1.0_prec,     &
                       mass = 1.0_prec,     &
                       !cubic_coupling = 390.0_prec, &
                       cubic_coupling = 5274.0_prec, &
                       potential = harmonic_trap_2D,   &
                       boundary_conditions = periodic)

    psi = wf_schroedinger_real_2D(method)

  !  method_h = schroedinger_hermite_real_2D(nx = 100, omega_x = 1.0_prec,&
  !                                          ny = 100, omega_y = 1.0_prec,&
! !                                          cubic_coupling = 390.0_prec)
   !                                         cubic_coupling = 5274.0_prec)
 !   psi_h = wf_schroedinger_hermite_real_2D(method_h)

    

 !   psi%u = 1.0_prec
 !   call psi%print_local_imaginary_time_orders(dt0 = 0.1_prec,     &
 !                                rows = 15, &
  !                               start_with_B=.true., &
  !                               extrapolation_order = 2)


    call psi%compute_groundstate(dt0 = 0.025_prec,     &
                                 tol = 1e-8_prec,    &
                                 max_iters = 10000,  &
                                 start_with_B=.true., &
                                 extrapolation_order = 2)
    call psi%save("xxx.h5")


!    psi_h%u = 1.0_prec
!    call psi_h%print_local_imaginary_time_orders(dt0 = 0.1_prec,     &
!                                 rows = 15, &
!                                 start_with_B=.true., &
!                                 extrapolation_order = 2)
!
!
!    call psi_h%compute_groundstate(dt0 = 0.025_prec,     &
!                                 tol = 1e-6_prec,    &
!                                 max_iters = 10000,  &
!                                 start_with_B=.true., &
!                                 extrapolation_order = 2)
!    call psi_h%save("xxx_h.h5")



    call finalize_tssm

#else
    type(schroedinger_1D) :: method 
    type(wf_schroedinger_1D) :: psi, psi0
    
    complex(prec) :: omega1, omega2
    complex(prec) :: scheme_A44c(7)
!    omega1 =  1.351207191959657634047688_prec ! 1/(2-2^(1/3))
!    omega2 = -1.702414383919315268095376_prec 
    omega1 = 1.0_prec/(2.0_prec - 2.0_prec**(1.0_prec/3.0_prec)*exp(cmplx(0.0_prec, 2.0_prec*pi/3.0_prec, prec)))
    omega2 = 1.0_prec -2.0_prec*omega1
    scheme_A44c = (/ 0.5_prec*omega1, omega1, 0.5_prec*(omega1+omega2), omega2, 0.5_prec*(omega1+omega2), omega1, 0.5_prec*omega1 /)  



    call initialize_tssm


    method = schroedinger_1D(nx = 512,     &
                       xmin = -20.0_prec,  &
                       xmax =  20.0_prec,  &
!                       ny = 512, &
!                       ymin = -20.0_prec,  &
!                       ymax =  20.0_prec,  &
                       hbar = 1.0_prec,     &
                       mass = 1.0_prec,     &
                       cubic_coupling = 390.0_prec, &
                       !cubic_coupling = 5274.0_prec, &
                       potential = harmonic_trap_1D,   &
                       boundary_conditions = dirichlet)
    
    psi =  wf_schroedinger_1D(method)

!    psi0 =  wf_schroedinger_1D(method)
    !psi%u = 1.0_prec
    call psi%rset(init)
!    call psi0%rset(init)
!    print *, "N", psi%norm2()
!    call psi%selfconsistent_nonlinear_step((0.001_prec,0.001_prec), (0.0005_prec,.0005_prec))
!    call psi%selfconsistent_nonlinear_step((-0.001_prec,-0.001_prec), (-0.0005_prec,-.0005_prec))
!    print *, "D", psi%distance(psi0)
!    call psi%save("xxx.h5")



    call psi%print_local_imaginary_time_orders(dt0 = 0.1_prec,     &
                                 rows = 15, &
                                 start_with_B=.false., &
                                 splitting_scheme = scheme_A44c, &
                                 extrapolation_order = 2)


!    call psi%compute_groundstate(dt0 = 0.05_prec,     &
!                                 tol = 1e-10_prec,    &
!                                 max_iters = 10000,  &
!                                 start_with_B=.true., &
!                                 splitting_scheme = scheme_A44c, &
!                                 extrapolation_order = 1)

!    call psi%save("xxx.h5")


    call finalize_tssm


#endif 

end program test_gross_pitaevskii
