program test_compare_space_discretizations
    use tssm
    implicit none

    integer, parameter :: ex0 = 10 
    integer, parameter :: ex1 = 17 
    integer :: ex

    real(prec), parameter :: g0 = 390.3654872824144

    real(prec), parameter :: xmin = -1000.0_prec
    real(prec), parameter :: xmax = +1000.0_prec
    integer :: nx 

    !smaller domain for groundstate:
    real(prec), parameter :: xmin_gs = xmin/8.0_prec
    real(prec), parameter :: xmax_gs = xmax/8.0_prec
    integer :: nx_gs

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
    real(prec) :: ell ! for periodic potential
    
    character(len=2) :: s
    character(len=20) :: filename_groundstate
    character(len=20) :: filename_V_disorder

    call initialize_tssm

#if 1    

    do ex = ex0, ex1
         nx = 2**ex
         print *, "*** ex =", ex, "  nx =", nx
         nx_gs = nx/8
         write (s,'(I2)') ex
         filename_groundstate = "groundstate_e" // s // ".h5"
         filename_V_disorder = "V_disorder_e" // s // ".h5"

         call compute_groundstate
    
         V0 = 0.2_prec*E
         !xi = 1/sqrt(4.0_prec*mu)
         xi = 1.0_prec/sqrt(8.0_prec*E)   
         kmax = 1.0_prec/xi  
         sigma = 0.7_prec*xi
         ell = 5.8_prec * xi

         if (this_proc==0) then
             print *, "mu=", mu
             print *, "E=", E
             print *, "xi=", xi
             print *, "kmax=", kmax
             print *, "sigma=", sigma 
             print *, "V0=", V0
             print *, "ell=", ell
         endif

         !call compute_disorder_potential
         
    end do
#endif 

    call propagate_gpe_parallel_and_compare

    call finalize_tssm

contains

subroutine propagate_gpe_parallel_and_compare
    use tssm_schroedinger
    use tssm_splitting_schemes
    type(schroedinger_1D) :: method(ex0:ex1) 
    type(wf_schroedinger_1D) :: psi(ex0:ex1)

    character(len=2) :: s
    character(len=30) :: filename

    integer(kind=8), save :: calc 
    real :: calc_time         


    real(prec) :: dt, t
    integer :: step, steps
    real(prec) :: d(ex0:ex1-1)

    if (this_proc==0) then
       print *, "*** Propagate GPE and compare"
    endif

    call tick(calc)

    do ex = ex0, ex1
        nx = 2**ex
        nx_gs = nx/8
        write (s,'(I2)') ex
        method(ex) = schroedinger_1D(nx = nx, &
                   xmin = xmin,       &
                   xmax = xmax,       &
                   hbar = 1.0_prec,          &
                   mass = 1.0_prec,          &
                   cubic_coupling = g0, &
                   boundary_conditions = periodic)
#if 0                   
        filename = "V_disorder_e" // s // ".h5"
        call method(ex)%load_potential(trim(filename))
#else
        call method(ex)%set_potential(periodic_potential_1D)
#endif

        psi(ex) =  wf_schroedinger_1D(method(ex))     
        filename = "groundstate_e" // s // ".h5"
        call psi(ex)%load(trim(filename))  
    end do

    dt = 2.0e-3_prec
    steps = 12000
    t = 0.0_prec
    step = 0

        call psi(ex1)%to_real_space
        do ex = ex0, ex1-1
            call psi(ex)%to_real_space
            psi(ex)%u = psi(ex)%u-psi(ex1)%u(lbound(psi(ex1)%u,1)+2**(ex1-ex)-1::2**(ex1-ex))
            d(ex) = psi(ex)%norm2()
            psi(ex)%u = psi(ex)%u+psi(ex1)%u(lbound(psi(ex1)%u,1)+2**(ex1-ex)-1::2**(ex1-ex))
        end do
        calc_time = tock(calc)
        write (*, '(I5,F10.2,20E23.15)') step, calc_time, t, (d(ex), ex=ex0,ex-1)

    do step = 1, steps
    
        do ex = ex0, ex1
            call psi(ex)%step(dt, t, & 
                          splitting_scheme = coeffs_15A, &
                          operator_sequence = "AB")
        end do          


        call psi(ex1)%to_real_space
        do ex = ex0, ex1-1
            call psi(ex)%to_real_space
            psi(ex)%u = psi(ex)%u-psi(ex1)%u(lbound(psi(ex1)%u,1)+2**(ex1-ex)-1::2**(ex1-ex))
            d(ex) = psi(ex)%norm2()
            psi(ex)%u = psi(ex)%u+psi(ex1)%u(lbound(psi(ex1)%u,1)+2**(ex1-ex)-1::2**(ex1-ex))
        end do
        t = t + dt
        calc_time = tock(calc)
        write (*, '(I5,F10.2,20E23.15)') step, calc_time, t, (d(ex), ex=ex0,ex-1)
                          
        if (step==6000) then
            do ex = ex0, ex1
                write (s,'(I2)') ex
                !filename = "solution_t12_e" // s // ".h5"
                filename = "solution_t12_periodic_e" // s // ".h5"
                call psi(ex)%save(trim(filename))
            end do    
        else if (step==12000) then
            do ex = ex0, ex1
                write (s,'(I2)') ex
                !filename = "solution_t24_e" // s // ".h5"
                filename = "solution_t24_periodic_e" // s // ".h5"
                call psi(ex)%save(trim(filename))
            end do    
        end if
   end do     

end subroutine propagate_gpe_parallel_and_compare


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

#if 1    
    call psi%compute_groundstate(dt0 = 0.05_prec, &
                   tol = 1e-10_prec,    &
                   max_iters = 10000,   &
                   extrapolation_order = 2)

    call psi%save(trim(filename_groundstate))                       
#else    
    call psi%load(trim(filename_groundstate))                       
#endif    

    !Compute chemical potential:
    !mu = psi%kinetic_energy() + psi%potential_energy() + 2.0_prec*psi%interaction_energy() 
    mu = psi%kinetic_energy()  + 2.0_prec*psi%interaction_energy() 
    ! ? mit Falle besser mit potential energy to be included, für  Rechnungen nachvollziehen besser da 
   
    !Compute energy:
    E = psi%kinetic_energy() + psi%interaction_energy()
    !Note: energy of groundstate (computed with harmonic pot) but potential
    !energy not included.

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

    call m%save_potential(trim(filename_V_disorder))

end subroutine compute_disorder_potential


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
        V = 0.5_prec * x**2
    end function harmonic_trap_1D

    function periodic_potential_1D(x) result(V)
        real(kind=prec), intent(in) :: x
        real(kind=prec) :: V
        real(kind=prec), parameter :: pi = 3.1415926535897932384626433832795028841971693993751_prec
        V = V0 * cos(2.0_prec*pi/ell * x)
    end function periodic_potential_1D

   

end program test_compare_space_discretizations
