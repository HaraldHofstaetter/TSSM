#define _ADAPTIVE_
program test_paper_fig_8
    use tssm
    implicit none

    real(prec), parameter :: g0 = 390.3654872824144

    real(prec), parameter :: xmin = -1000.0_prec
    real(prec), parameter :: xmax = +1000.0_prec
    integer, parameter :: nx = 2**17

    !smaller domain for groundstate:
    real(prec), parameter :: xmin_gs = xmin/8.0_prec
    real(prec), parameter :: xmax_gs = xmax/8.0_prec
    integer, parameter :: nx_gs = nx/8

    !parameters for disorder potential:
!    real(prec), parameter :: x_first_gaussian = xmin
!    real(prec), parameter :: x_last_gaussian = xmax
    real(prec), parameter :: x_first_gaussian = -999.99060221295724_prec
    real(prec), parameter :: x_last_gaussian = 999.99060221338993_prec
    integer, parameter :: number_of_gaussians = 20000
    real(prec) :: V0     ! to be computed
    real(prec) :: sigma  ! to be computed

    !perturbance parameter for intitial states
    real(prec), parameter :: alpha = 1.0E-4_prec 

    real(prec) :: mu, xi, kmax, E 

    call initialize_tssm

#if 1    

#if 0
    mu=   27.975484637775107     
    E=   13.996353367380795     
    xi=   9.4503426885274647E-002
    kmax=   10.581626856917907     
    sigma=   6.6152398819692246E-002
    V0=   2.7992706734761592     
#else

    call compute_groundstate
    
    V0 = 0.2_prec*E
    !xi = 1/sqrt(4.0_prec*mu)
    xi = 1.0_prec/sqrt(8.0_prec*E)   
    kmax = 1.0_prec/xi  
    sigma = 0.7_prec*xi

    if (this_proc==0) then
        print *, "mu=", mu
        print *, "E=", E
        print *, "xi=", xi
        print *, "kmax=", kmax
        print *, "sigma=", sigma 
        print *, "V0=", V0
    endif
#endif


    call compute_disorder_potential
#endif    

    call propagate_gpe

    call finalize_tssm

contains   

subroutine compute_groundstate
    use tssm_schroedinger
    type(schroedinger_real_1D) :: method 
    type(wf_schroedinger_real_1D) :: psi

    integer :: ix
    real(prec) :: x

    if (this_proc==0) then
       print *, "*** Computing groundstate for BEC..."
    endif

    method = schroedinger_real_1D(nx = nx_gs, &
                   xmin = xmin_gs,       &
                   xmax = xmax_gs,       &
                   hbar = 1.0_prec,          &
                   mass = 1.0_prec,          &
                   cubic_coupling = g0, &
                   potential = harmonic_trap_1D,    &
                   boundary_conditions = periodic)

    psi =  wf_schroedinger_real_1D(method)     

    call psi%compute_groundstate(dt0 = 0.05_prec, &
                   tol = 1e-10_prec,    &
                   max_iters = 10000,   &
                   extrapolation_order = 2)

    call psi%save("groundstate.h5")                       

    !Compute chemical potential:
    !mu = psi%kinetic_energy() + psi%potential_energy() + 2.0_prec*psi%interaction_energy() 
    mu = psi%kinetic_energy()  + 2.0_prec*psi%interaction_energy() 
    ! ? mit Falle besser mit potential energy to be included, für  Rechnungen nachvollziehen besser da 
   
    !Compute energy:
    E = psi%kinetic_energy() + psi%interaction_energy()
    !Note: energy of groundstate (computed with harmonic pot) but potential
    !energy not included.

#if 1
    !Compute perturbed initial states
    call psi%to_real_space
    do ix=method%g%n1min, method%g%n1max
       x = method%g%xmin + ix*method%g%dx
       psi%u(ix) = psi%u(ix) * (1.0_prec-alpha*x)
    end do   
    call psi%normalize
    call psi%save("groundstate_perturbed_1.h5")

    call psi%load("groundstate.h5")
    call psi%to_real_space
    do ix=method%g%n1min, method%g%n1max
       x = method%g%xmin + ix*method%g%dx
       psi%u(ix) = psi%u(ix) * (1.0_prec+alpha*x)
    end do   
    call psi%normalize    
    call psi%save("groundstate_perturbed_2.h5")
#endif    

 !p [-10:+10] "< h5totxt -d psi_real groundstate_perturbed_1.h5" u (-125+$0*250/2**12):($1) w l,  "< h5totxt -d psi_real groundstate_perturbed_2.h5" u (-125+$0*250/2**12):($1) w l
      
    call psi%finalize
       
end subroutine compute_groundstate


subroutine compute_disorder_potential    
    use tssm_schroedinger
    use tssm_disorder_potential
    implicit none
    type(schroedinger_1D) :: m
    real(prec) :: mean, dev

    if (this_proc==0) then
       print *, "*** Computing disordered potential..."
    endif

    m = schroedinger_1D(nx = nx, xmin = xmin, xmax = xmax, & 
             potential = harmonic_trap_1D) ! dummy potential to be replaced by disorder potential

    call compute_disorder_potential_ala_iva_1D(m%g, m%V, &
                   x_first_gaussian = x_first_gaussian, &
                   x_last_gaussian = x_last_gaussian, &
                   number_of_gaussians = number_of_gaussians, &
                   sigma = sigma, &
                   f = 10)        

    call m%save_potential("V_disorder_unnormalized.h5")
    call get_mean_value_and_standard_deviation_1D(m%g, m%V, mean, dev)
    if (this_proc==0) then
        print *, "V_disorder (unnormalized): mean, dev =", mean, dev
    end if    

    !normalize potential
    m%V = (V0/dev)*(m%V - mean)

    call get_mean_value_and_standard_deviation_1D(m%g, m%V, mean, dev)
    if (this_proc==0) then
        print *, "V_disorder (normalized): mean, dev =", mean, dev
    end if    

    call m%save_potential("V_disorder.h5")

! p [999:1000] "< h5totxt /home/hofi/TSSM/OPENMP/V_disorder.h5" u (-1000+($0+1)*2000/2**16):($1) w lp, "potential_wi.dat" u 1:3 w lp

end subroutine compute_disorder_potential


subroutine propagate_gpe
    use tssm_schroedinger
    use tssm_splitting_schemes
    type(schroedinger_1D) :: method 
    type(wf_schroedinger_1D) :: psi

    real(prec) :: dt, tol

    if (this_proc==0) then
       print *, "*** Propagate GPE ... disturbation"
    endif


    method = schroedinger_1D(nx = nx, &
                   xmin = xmin,       &
                   xmax = xmax,       &
                   hbar = 1.0_prec,          &
                   mass = 1.0_prec,          &
                   cubic_coupling = g0, &
                   boundary_conditions = periodic)
   
   call method%load_potential("V_disorder.h5")

   psi =  wf_schroedinger_1D(method)     

   !call psi%load("groundstate.h5")  
   call psi%load("groundstate_perturbed_1.h5")  

   !tol = 1e-9_prec
   tol = 1e-8_prec

#ifdef _ADAPTIVE_
   dt = 2.0e-3_prec
   call psi%run_adaptive(dt, (/ 0.0_prec, 12.0_prec /), &
                tol = tol, &
                splitting_scheme = coeffs_15A, &
                controller_scheme = coeffs_15, &
                order = 4, &
                operator_sequence = "AB", &
                solution_out = solution_out_schroedinger_1D)

   !call psi%save("solution_t12_pert1.h5")    
   call psi%save("solution_t12_pert1_tol1e-8.h5")    

   call psi%run_adaptive(dt, (/ 12.0_prec, 24.0_prec /), &
                tol = tol, &
                splitting_scheme = coeffs_15A, &
                controller_scheme = coeffs_15, &
                order = 4, &
                operator_sequence = "AB", &
                solution_out = solution_out_schroedinger_1D)

   !call psi%save("solution_t24_pert1.h5")    
   call psi%save("solution_t24_pert1_tol1e-8.h5")    

   call psi%run_adaptive(dt, (/ 24.0_prec, 0.0_prec /), &
                tol = tol, &
                splitting_scheme = coeffs_15A, &
                controller_scheme = coeffs_15, &
                order = 4, &
                operator_sequence = "AB", &
                solution_out = solution_out_schroedinger_1D)

   !call psi%save("solution_backward_t24_to_t0_pert1.h5")  
   call psi%save("solution_backward_t24_to_t0_pert1_tol1e-8.h5")  

#else

   dt = 2.0e-3_prec
   call psi%run(dt, (/ 0.0_prec, 12.0_prec /), &
                splitting_scheme = coeffs_15A, &
                operator_sequence = "AB", &
                solution_out = solution_out_schroedinger_1D)

   call psi%save("solution_t12_csz.h5")    

   call psi%run(dt, (/ 12.0_prec, 24.0_prec /), &
                splitting_scheme = coeffs_15A, &
                operator_sequence = "AB", &
                solution_out = solution_out_schroedinger_1D)

   call psi%save("solution_t24_csz.h5")    

   call psi%run(dt, (/ 24.0_prec, 0.0_prec /), &
                splitting_scheme = coeffs_15A, &
                operator_sequence = "AB", &
                solution_out = solution_out_schroedinger_1D)

   call psi%save("solution_backward_t24_to_t0_csz.h5")  

#endif

   call psi%finalize
        
end subroutine propagate_gpe




    subroutine compute_disorder_potential_ala_iva_1D(g, V, x_first_gaussian, x_last_gaussian, &
                                                     number_of_gaussians, sigma, f) 
        use tssm_schroedinger
        use tssm_disorder_potential
        class(grid_equidistant_1D) :: g
        real(kind=prec), intent(inout) :: V(g%n1min:g%n1max)
        real(kind=prec), intent(in) :: x_first_gaussian
        real(kind=prec), intent(in) :: x_last_gaussian 
        integer, intent(in) :: number_of_gaussians 
        real(kind=prec), intent(in) :: sigma
        integer, intent(in) :: f

        INTEGER, PARAMETER :: I4B=SELECTED_INT_KIND(9)
        integer(kind=I4B) :: idum

        integer :: shell, shells
        integer :: ix
        real(prec) :: A, x_center, dx

        idum= -4309456 ! eine Hausnummer
        dx = (x_last_gaussian-x_first_gaussian)/real(number_of_gaussians-1, kind=prec)

        V = 0.0_prec

        do ix=1, number_of_gaussians
            x_center = x_first_gaussian + (ix-1)*dx
            A = ran(idum)
!           call add_gaussian_1D(g, V, A-0.5_prec, x_center, sigma, f)
            call add_gaussian_1D(g, V, A, x_center, sigma, f)
        end do

    end subroutine compute_disorder_potential_ala_iva_1D


    FUNCTION ran(idum)
    IMPLICIT NONE
    INTEGER, PARAMETER :: I4B=SELECTED_INT_KIND(9)
    INTEGER(I4B), INTENT(INOUT) :: idum
    REAL :: ran
       !Minimal random number generator of Park and Miller combined with a Marsaglia shift
       !sequence. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
       !values). This fully portable, scalar generator has the N4traditionalN! (not Fortran 90) calling
       !sequence with a random deviate as the returned function value: call with idum a negative
       !integer to initialize; thereafter, do not alter idum except to reinitialize. The period of this
       !generator is about 3.1 NW 1018 .
    INTEGER(I4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
    REAL, SAVE :: am
    INTEGER(I4B), SAVE :: ix=-1,iy=-1,k
    
    !Initialize.
    IF (idum <= 0 .OR. iy < 0) THEN
        am=NEAREST(1.0,-1.0)/IM
        iy=IOR(IEOR(888889999,ABS(idum)),1)
        ix=IEOR(777755555,ABS(idum))
        !Set idum positive.
        idum=ABS(idum)+1
    END IF
    !Marsaglia shift sequence with period 232 $(B!](B 1.
    ix=IEOR(ix,ISHFT(ix,13))
    ix=IEOR(ix,ISHFT(ix,-17))
    ix=IEOR(ix,ISHFT(ix,5))
    !Park-Miller sequence by SchrageN"s method, period 231 $(B!](B 2.
    k=iy/IQ
    iy=IA*(iy-k*IQ)-IR*k
    IF (iy < 0) iy=iy+IM
    !Combine the two generators with masking to ensure nonzero value.
    ran=am*IOR(IAND(IM,IEOR(ix,iy)),1)
    
    END FUNCTION ran



    function harmonic_trap_1D(x) result(V)
        real(kind=prec), intent(in) :: x
        real(kind=prec) :: V
        !V = 0.5_prec * omega**2 * x**2
        V = 0.5_prec * x**2
    end function harmonic_trap_1D

end program test_paper_fig_8
