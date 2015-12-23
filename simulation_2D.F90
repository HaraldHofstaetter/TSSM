program simulation_2D    
    use tssm_schroedinger
    implicit none

    real(prec), parameter :: N =2.3e5_prec
    real(prec), parameter :: cubic_coupling = 5274.41_prec !System 3
    real(prec), parameter :: cubic_coupling_dummy = 2.0_prec * cubic_coupling

    real(prec), parameter :: xmin_gs = -25.0_prec
    real(prec), parameter :: xmax_gs = +25.0_prec
    integer, parameter :: nx_gs = 512
    real(prec), parameter :: ymin_gs = -25.0_prec
    real(prec), parameter :: ymax_gs = +25.0_prec
    integer, parameter :: ny_gs = 512

    real(prec), parameter :: xmin = 4*xmin_gs
    real(prec), parameter :: xmax = 4*xmax_gs
    integer, parameter :: nx = 4*nx_gs 
    real(prec), parameter :: ymin = 4*ymin_gs
    real(prec), parameter :: ymax = 4*ymax_gs
    integer, parameter :: ny = 4*ny_gs 

    real(prec) :: mu, xi, kmax, sigma, V0, E_gpe_tot

    call initialize_tssm

    call compute_groundstate_gpe

    xi = 1/sqrt(4.0_prec*mu)
    ! kmax = 1.0_prec/xi  ! for 1D
    kmax = sqrt(0.5_prec)/xi ! for 2D
    sigma = 0.5_prec/kmax
    if (this_proc==0) then
        print *, "mu", mu
        print *, "xi", xi
        print *, "kmax", kmax
        print *, "sigma", sigma 
        print *, "V0", V0
    endif

    call compute_disordered_potential

    call compute_groundstate_dummy

    call compute_initial_state_for_SE

    call propagate_SE

    call finalize_tssm

contains     

  subroutine compute_groundstate_gpe
        use tssm_schroedinger
        type(schroedinger_real_2D) :: method 
        type(wf_schroedinger_real_2D) :: psi

        if (this_proc==0) then
           print *, "*** Computing groundstate for BEC..."
        endif

        method = schroedinger_real_2D(nx = nx_gs, &
                       xmin = xmin_gs,       &
                       xmax = xmax_gs,       &
                       ny = ny_gs,           &
                       ymin = ymin_gs,       &
                       ymax = ymax_gs,       &
                       hbar = 1.0_prec,          &
                       mass = 1.0_prec,          &
                       cubic_coupling = cubic_coupling, &
                       potential = harmonic_trap_2D,    &
                       boundary_conditions = dirichlet)

       psi =  wf_schroedinger_real_2D(method)     

       call psi%compute_groundstate(dt0 = 0.05_prec, &
                       tol = 1e-10_prec,    &
                       max_iters = 10000,   &
                       extrapolation_order = 2)

       call psi%save("groundstate_gpe.h5")                       

       !Chemical potential:

       mu = psi%kinetic_energy() + psi%potential_energy() + 2.0_prec*psi%interaction_energy()
       V0 = 0.2_prec*(psi%kinetic_energy() + psi%interaction_energy())
      
       call psi%finalize
       
    end subroutine compute_groundstate_gpe


    subroutine compute_disordered_potential 
        use tssm_schroedinger
        use tssm_disorder_potential


        type(schroedinger_2D) :: m
        type(wf_schroedinger_2D) :: psi
        real(prec) :: mean, dev, mu, xi

        if (this_proc==0) then
           print *, "*** Computing disordered potential..."
        endif

         m = schroedinger_2D(nx = nx,     &
                       xmin = xmin,       &
                       xmax = xmax,       &
                       ny = ny,           &
                       ymin = ymin,       &
                       ymax = ymax,       &
                       hbar = 1.0_prec,          &
                       mass = 1.0_prec,          &
                       cubic_coupling = cubic_coupling, &
                       potential = harmonic_trap_2D,    &
                       boundary_conditions = periodic)

        call compute_disorder_potential_2D(m%g, m%V, &
                   dx = 0.1_prec, dy = 0.1_prec, &
                   sigma = sigma, &
                   f = 10)        


       call m%save_potential("V_disorder_unnormalized.h5")
       call get_mean_value_and_standard_deviation_2D(m%g, m%V, mean, dev)
       if (this_proc==0) then
           print *, "mean, dev =", mean, dev
       end if    

       !normalize potential
       m%V = (V0/dev)*(m%V - mean)

       call get_mean_value_and_standard_deviation_2D(m%g, m%V, mean, dev)
       if (this_proc==0) then
           print *, "mean, dev =", mean, dev
       end if    

       call m%save_potential("V_disorder.h5")

       if (this_proc==0) then
           print *, "*** Computing total energy for BEC ..."
       endif

       psi =  wf_schroedinger_2D(m)     
       call psi%load("groundstate_gpe.h5") 
       E_gpe_tot = psi%kinetic_energy() + psi%potential_energy() + psi%interaction_energy()

       if (this_proc==0) then
           print  * ,"E_gpe_tot", E_gpe_tot
       endif     

       call psi%finalize

    end subroutine compute_disordered_potential

    subroutine compute_groundstate_dummy
        use tssm_schroedinger
        type(schroedinger_real_2D) :: method 
        type(wf_schroedinger_real_2D) :: psi

        if (this_proc==0) then
           print *, "*** Computing dummy groundstate"
        endif

        method = schroedinger_real_2D(nx = nx_gs, &
                       xmin = xmin_gs,       &
                       xmax = xmax_gs,       &
                       ny = ny_gs,           &
                       ymin = ymin_gs,       &
                       ymax = ymax_gs,       &
                       hbar = 1.0_prec,          &
                       mass = 1.0_prec,          &
                       cubic_coupling = cubic_coupling_dummy, & !!!<<<------!!!
                       potential = harmonic_trap_2D,    &
                       boundary_conditions = dirichlet)

       psi =  wf_schroedinger_real_2D(method)     

       call psi%compute_groundstate(dt0 = 0.05_prec, &
                       tol = 1e-10_prec,    &
                       max_iters = 10000,   &
                       extrapolation_order = 2)

       call psi%save("groundstate_dummy.h5")                       

       call psi%finalize
       
    end subroutine compute_groundstate_dummy



    subroutine compute_initial_state_for_SE
        use tssm_schroedinger
        use tssm_splitting_schemes
        type(schroedinger_2D) :: method 
        type(wf_schroedinger_2D) :: psi
        real(prec) :: dt

        if (this_proc==0) then
           print *, "*** Computing initial state for SE"
           print *, "*** propagate until E_kin>=E_gpe_tot =", E_gpe_tot
        endif

        method = schroedinger_2D(nx = nx,     &
                 xmin = xmin,       &
                 xmax = xmax,       &
                 ny = ny,           &
                 ymin = ymin,       &
                 ymax = ymax,       &
                 hbar = 1.0_prec,          &
                 mass = 1.0_prec,          &
                 cubic_coupling = cubic_coupling_dummy, & !!!<<<------!!
                 boundary_conditions = periodic)

         psi = wf_schroedinger_2D(method)
         
         call psi%load("groundstate_dummy.h5") 
         dt = 1.0e-3_prec
         call psi%run_adaptive(dt = dt, &
                    t0=0.0_prec, &
                    tend = 1e5_prec, &   ! dummy value, 
                    tol = 1e-7_prec, &
                    splitting_scheme = coeffs_15, &
                    associated_scheme = coeffs_15A, &
                    order = 4, &
                    solution_out = solution_out_for_initial_state)

         call psi%save("initial_state_se.h5")              

        call psi%finalize
    end subroutine compute_initial_state_for_SE



     subroutine solution_out_for_initial_state(psi, t, dt, flag, abort_flag)
        use tssm_schroedinger
        class(wave_function), intent(inout) :: psi
        real(kind=prec), intent(in) :: t 
        real(kind=prec), intent(in) :: dt 
        character, intent(in) :: flag
        logical, intent(out) :: abort_flag
        integer, save :: step
        real(kind=prec) :: E_kin, E_int, E_tot
        real :: current_time
        real, save :: start_time

        select type (psi); class is (wf_schroedinger_2D)

        if (flag=="i" .or. flag=="I") then ! init =======================
            call cpu_time(start_time) 
            step = 0
            if (this_proc==0) then
              write (*,*), "                       t                  dt               E_kin" &
                           //"               E_int               E_tot"
              write (*,*), "------------------------------------------------------------------------------------"
            end if           
            abort_flag = .false.
        end if   

        if (flag=="e" .or. flag=="E") then ! execute =====
            step = step + 1
            E_kin = psi%kinetic_energy()
            E_int = psi%interaction_energy()
            E_tot = E_kin + E_int
            call cpu_time(current_time) 
            if (this_proc==0) then
               write (*, '(I4,5E20.12,F10.2)') step, t, dt,  E_kin, E_int, E_tot, current_time-start_time
            end if   
            abort_flag = E_kin>=E_gpe_tot 
        end if
    
        if (flag=="f" .or. flag=="F") then ! finish ====================
            abort_flag = .false.
        end if   

        end select
    end subroutine solution_out_for_initial_state
 

     subroutine propagate_SE
        use tssm_schroedinger
        use tssm_splitting_schemes
        type(schroedinger_2D) :: method 
        type(wf_schroedinger_2D) :: psi
        real(prec) :: dt
        real(prec) :: times(51) 
        integer :: j

        if (this_proc==0) then
           print *, "*** Propagate SE"
        endif

        method = schroedinger_2D(nx = nx,     &
                 xmin = xmin,       &
                 xmax = xmax,       &
                 ny = ny,           &
                 ymin = ymin,       &
                 ymax = ymax,       &
                 hbar = 1.0_prec,          &
                 mass = 1.0_prec,          &
                 cubic_coupling = 0.0_prec, & !!!<<<------!!
                 boundary_conditions = periodic)

         call method%load_potential("V_disorder.h5")

         psi = wf_schroedinger_2D(method)
         
         call psi%load("initial_state_se.h5")              
         dt = 1.0e-3_prec

!         call psi%run_adaptive(dt = dt, &
!                    t0 = 0.0_prec, &
!                    tend = 500.0_prec, &   ! dummy value, 
!                    tol = 1e-7_prec, &
!                    splitting_scheme = coeffs_15, &
!                    associated_scheme = coeffs_15A, &
!                    order = 4, &
!                    solution_out = solution_out)
!
!         call psi%save("solution_se_t500.h5")              

         times = (/ (10.0_prec*real(j, prec), j=0,50) /)

         call psi%run_adaptive_2(dt = dt, &
                    times = times, &
                    tol = 1e-7_prec, &
                    splitting_scheme = coeffs_15, &
                    associated_scheme = coeffs_15A, &
                    order = 4, &
                    solution_out = solution_out, &
                    solution_out_2 = solution_out_2)

        call psi%finalize
    end subroutine propagate_SE



     subroutine solution_out(psi, t, dt, flag, abort_flag)
        use tssm_schroedinger
        class(wave_function), intent(inout) :: psi
        real(kind=prec), intent(in) :: t 
        real(kind=prec), intent(in) :: dt 
        character, intent(in) :: flag
        logical, intent(out) :: abort_flag
        integer, save :: step
        real(kind=prec) :: E_kin, E_pot, E_int, E_tot, x_mean, x_dev, y_mean, y_dev
        real :: current_time
        real, save :: start_time = -1.0

        select type (psi); class is (wf_schroedinger_2D)

        if (flag=="i" .or. flag=="I") then ! init =======================
            if (start_time<0.0) then
                call cpu_time(start_time) 
            end if    
            step = 0
            if (this_proc==0) then
              write (*,*), "                       t                  dt               E_kin" &
                       //"               E_pot               E_int             E_tot          x_mean      x_dev elapsed_time"
              write (*,*), "------------------------------------------------------------------------------------"
            end if           
        end if   

        if (flag=="e" .or. flag=="E") then ! execute =====
            step = step + 1
            E_kin = psi%kinetic_energy()
            call psi%get_realspace_observables(E_pot, E_int, x_mean, x_dev, y_mean, y_dev)

            E_int = psi%interaction_energy()
            E_tot = E_kin + E_pot + E_int
            call cpu_time(current_time) 
            if (this_proc==0) then
               write (*, '(I6,11E20.12,F10.2)') step, t, dt, E_kin, E_pot, E_int, E_tot, &
                      x_mean, x_dev, y_mean, y_dev, sqrt(x_dev**2+y_dev**2), current_time-start_time
            end if   
        end if
    
        if (flag=="f" .or. flag=="F") then ! finish ====================
        end if   

        end select
        abort_flag = .false.
    end subroutine solution_out

    subroutine solution_out_2(psi, t, dt, flag, abort_flag)
        use tssm_schroedinger
        use tssm_hdf5
        class(wave_function), intent(inout) :: psi
        real(kind=prec), intent(in) :: t 
        real(kind=prec), intent(in) :: dt 
        character, intent(in) :: flag
        logical, intent(out) :: abort_flag
        integer, save :: step 
        character(len=4) :: s


        if (flag=="i" .or. flag=="I") then ! init =======================
            step = 0
        end if   

        if (flag=="e" .or. flag=="E") then ! execute =====
            step = step + 1
            write (s,'(I4.4)') step
            call psi%save("sol" // s // ".h5")
            call hdf5_write_attributes("sol" // s // ".h5" , dnames = (/ "t" /), dvalues = (/ t /) )
        end if   

        if (flag=="f" .or. flag=="F") then ! finish ====================
        end if   

        abort_flag = .false.
    end subroutine solution_out_2



  function harmonic_trap_2D(x,y) result(V)
        real(kind=prec), intent(in) :: x,y
        real(kind=prec) :: V
        V = 0.5_prec * (x**2 + y**2)
  end function harmonic_trap_2D





end program simulation_2D
