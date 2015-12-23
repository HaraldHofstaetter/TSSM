program simulation_1D    
    use tssm_schroedinger
    implicit none

    real(prec), parameter :: N =1.e4_prec
    real(prec), parameter :: cubic_coupling = 406.89_prec !System 3
    real(prec), parameter :: cubic_coupling_dummy = 2.0_prec * cubic_coupling

    real(prec), parameter :: xmin_gs = -100.0_prec
    real(prec), parameter :: xmax_gs = +100.0_prec
    integer, parameter :: nx_gs = 4096

    !real(prec), parameter :: xmin = 128*xmin_gs
    !real(prec), parameter :: xmax = 128*xmax_gs
    !integer, parameter :: nx = 128*nx_gs 

    real(prec), parameter :: xmin = 16*xmin_gs
    real(prec), parameter :: xmax = 16*xmax_gs
    integer, parameter :: nx = 16*nx_gs 

    real(prec) :: mu, xi, kmax, sigma, V0, E_gpe_tot

    call initialize_tssm

    call compute_groundstate_gpe

    xi = 1/sqrt(4.0_prec*mu)
    kmax = 1.0_prec/xi  
    ! kmax = sqrt(0.5_prec)/xi !(for 2D)
    sigma = 0.5_prec/kmax
    if (this_proc==0) then
        print *, "mu", mu
        print *, "xi", xi
        print *, "kmax", kmax
        print *, "sigma", sigma 
        print *, "V0", V0
    endif

!    call compute_disordered_potential

!    call compute_groundstate_dummy

!    call compute_initial_state_for_SE

!    call propagate_SE

    call propagate_GPE_free

!    call propagate_GPE

    call finalize_tssm

contains     

  subroutine compute_groundstate_gpe
        type(schroedinger_real_1D) :: method 
        type(wf_schroedinger_real_1D) :: psi

        if (this_proc==0) then
           print *, "*** Computing groundstate for BEC..."
        endif

        method = schroedinger_real_1D(nx = nx_gs, &
                       xmin = xmin_gs,       &
                       xmax = xmax_gs,       &
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

       call psi%save("groundstate_gpe.h5")                       

       !Chemical potential:

       mu = psi%kinetic_energy() + psi%potential_energy() + 2.0_prec*psi%interaction_energy() 
       ! ? mit Falle besser mit potential energy to be included, für  Rechnungen nachvollziehen besser da 
       V0 = 0.2_prec*(psi%kinetic_energy() + psi%interaction_energy()) 
      !Note: energy of groundstate (computed with harmonic pot) but potential
      !energy not included.
      
       call psi%finalize
       
    end subroutine compute_groundstate_gpe


    subroutine compute_disordered_potential 
        use tssm_disorder_potential


        type(schroedinger_1D) :: m
        type(wf_schroedinger_1D) :: psi
        real(prec) :: mean, dev, mu, xi

        if (this_proc==0) then
           print *, "*** Computing disordered potential..."
        endif

         m = schroedinger_1D(nx = nx,     &
                       xmin = xmin,       &
                       xmax = xmax,       &
                       hbar = 1.0_prec,          &
                       mass = 1.0_prec,          &
                       cubic_coupling = cubic_coupling, &
                       potential = harmonic_trap_1D,    &
                       boundary_conditions = periodic)

        call compute_disorder_potential_1D(m%g, m%V, &
                   dx = 0.16_prec,      &
                   sigma = sigma, &
                   f = 10)        


       call m%save_potential("V_disorder_unnormalized.h5")
       call get_mean_value_and_standard_deviation_1D(m%g, m%V, mean, dev)
       if (this_proc==0) then
           print *, "mean, dev =", mean, dev
       end if    

       !normalize potential
       m%V = (V0/dev)*(m%V - mean)

       call get_mean_value_and_standard_deviation_1D(m%g, m%V, mean, dev)
       if (this_proc==0) then
           print *, "mean, dev =", mean, dev
       end if    

       call m%save_potential("V_disorder.h5")

       if (this_proc==0) then
           print *, "*** Computing total energy for BEC ..."
       endif

       psi =  wf_schroedinger_1D(m)     
       call psi%load("groundstate_gpe.h5") 
       E_gpe_tot = psi%kinetic_energy() + psi%potential_energy() + psi%interaction_energy()

       if (this_proc==0) then
           print  * ,"E_gpe_tot", E_gpe_tot
       endif     

       call psi%finalize

    end subroutine compute_disordered_potential

    subroutine compute_groundstate_dummy
        use tssm_schroedinger
        type(schroedinger_real_1D) :: method 
        type(wf_schroedinger_real_1D) :: psi

        if (this_proc==0) then
           print *, "*** Computing dummy groundstate"
        endif

        method = schroedinger_real_1D(nx = nx_gs, &
                       xmin = xmin_gs,       &
                       xmax = xmax_gs,       &
                       hbar = 1.0_prec,          &
                       mass = 1.0_prec,          &
                       cubic_coupling = cubic_coupling_dummy, & !!!<<<------!!!
                       potential = harmonic_trap_1D,    &
                       boundary_conditions = periodic)

       psi =  wf_schroedinger_real_1D(method)     

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
        type(schroedinger_1D) :: method 
        type(wf_schroedinger_1D) :: psi
        real(prec) :: dt

        if (this_proc==0) then
           print *, "*** Computing initial state for SE"
           print *, "*** propagate until E_kin>=E_gpe_tot =", E_gpe_tot
        endif

        method = schroedinger_1D(nx = nx,     &
                 xmin = xmin,       &
                 xmax = xmax,       &
                 hbar = 1.0_prec,          &
                 mass = 1.0_prec,          &
                 cubic_coupling = cubic_coupling_dummy, & !!!<<<------!!
                 boundary_conditions = periodic)

         psi = wf_schroedinger_1D(method)
         
         call psi%load("groundstate_dummy.h5") 
         dt = 1.0e-3_prec
         call psi%run_adaptive(dt = dt, &
                    times = (/ 0.0_prec, &
                               1e5_prec /), &   ! dummy value, 
                    tol = 1e-7_prec, &
                    splitting_scheme = coeffs_15, &
                    controller_scheme = coeffs_15A, &
                    order = 4, &
                    solution_out = solution_out_for_initial_state)

         call psi%save("initial_state_se.h5")              

        call psi%finalize
    end subroutine compute_initial_state_for_SE



     subroutine solution_out_for_initial_state(psi, t, dt, flag, abort_flag)
        use tssm_common
        use tssm_schroedinger
        class(wave_function), intent(inout) :: psi
        real(kind=prec), intent(in) :: t 
        real(kind=prec), intent(in) :: dt 
        character, intent(in) :: flag
        logical, intent(out) :: abort_flag
        integer, save :: step
        real(kind=prec) :: E_kin, E_int, E_tot
        integer(kind=8), save :: calc
        real :: calc_time 

        select type (psi); class is (wf_schroedinger_1D)

        if (flag=="i" .or. flag=="I") then ! init =======================
            call tick(calc)
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
            calc_time = tock(calc)
            if (this_proc==0) then
               write (*, '(I4,5E20.12,F10.2)') step, t, dt,  E_kin, E_int, E_tot, calc_time
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
        type(schroedinger_1D) :: method 
        type(wf_schroedinger_1D) :: psi
        real(prec) :: dt
        real(prec) :: times(0:100) 
        integer :: j

        if (this_proc==0) then
           print *, "*** Propagate SE"
        endif

        method = schroedinger_1D(nx = nx,     &
                 xmin = xmin,       &
                 xmax = xmax,       &
                 hbar = 1.0_prec,          &
                 mass = 1.0_prec,          &
                 cubic_coupling = 0.0_prec, & !!!<<<------!!
                 boundary_conditions = periodic)

         call method%load_potential("V_disorder.h5")
         psi = wf_schroedinger_1D(method)
         
         call psi%load("initial_state_se.h5")  
         dt = 1.0e-3_prec

         times = (/ (10.0_prec*real(j, prec), j=0,100) /)

         call psi%run_adaptive(dt = dt, &
                    times = times, &
                    tol = 1e-7_prec, &
                    splitting_scheme = coeffs_15, &
                    controller_scheme = coeffs_15A, &
                    order = 4, &
                    solution_out = solution_out_se)

        call psi%finalize
    end subroutine propagate_SE


     subroutine propagate_GPE
        use tssm_schroedinger
        use tssm_splitting_schemes
        type(schroedinger_1D) :: method
        type(wf_schroedinger_1D) :: psi
        real(prec) :: dt
        real(prec) :: times(0:100)
        integer :: j

        if (this_proc==0) then
           print *, "*** Propagate GPE"
        endif

        method = schroedinger_1D(nx = nx,     &
                 xmin = xmin,       &
                 xmax = xmax,       &
                 hbar = 1.0_prec,          &
                 mass = 1.0_prec,          &
                 cubic_coupling = cubic_coupling, & !!!<<<------!!
                 boundary_conditions = periodic)

         call method%load_potential("V_disorder.h5")
         psi = wf_schroedinger_1D(method)

         call psi%load("groundstate_gpe.h5") 

         dt = 1.0e-3_prec

         times = (/ (10.0_prec*real(j, prec), j=0,100) /)

         call psi%run_adaptive(dt = dt, &
                    times = times, &
                    tol = 1e-7_prec, &
                    splitting_scheme = coeffs_15, &
                    controller_scheme = coeffs_15A, &
                    order = 4, &
                    solution_out = solution_out_gpe)

        call psi%finalize
    end subroutine propagate_GPE


     subroutine propagate_GPE_free
        use tssm_schroedinger
        use tssm_splitting_schemes
        type(schroedinger_1D) :: method
        type(wf_schroedinger_1D) :: psi
        real(prec) :: dt
        real(prec) :: times(0:100)
        integer :: j

        if (this_proc==0) then
           print *, "*** Propagate GPE (free)"
        endif

        method = schroedinger_1D(nx = nx,     &
                 xmin = xmin,       &
                 xmax = xmax,       &
                 hbar = 1.0_prec,          &
                 mass = 1.0_prec,          &
                 cubic_coupling = cubic_coupling, & !!!<<<------!!
                 boundary_conditions = periodic)
         !Note: no potential

         psi = wf_schroedinger_1D(method)

         call psi%load("groundstate_gpe.h5") 

         dt = 1.0e-3_prec

         times = (/ (10.0_prec*real(j, prec), j=0,100) /)

         !call psi%run_adaptive(dt = dt, &
         call psi%run(dt = dt, &
                    times = times, &
                    !tol = 1e-8_prec, &
                    splitting_scheme = coeffs_15, &
                    !controller_scheme = coeffs_15A, &
                    !order = 4, &
                    solution_out = solution_out_gpe_free)

        call psi%finalize
    end subroutine propagate_GPE_free


    subroutine solution_out_se(psi, t, dt, flag, abort_flag)
        use tssm_hdf5
        class(wave_function), intent(inout) :: psi
        real(kind=prec), intent(in) :: t 
        real(kind=prec), intent(in) :: dt 
        character, intent(in) :: flag
        logical, intent(out) :: abort_flag
        integer, save :: big_step 
        character(len=4) :: s

        if (flag=="i" .or. flag=="I") then ! init =====
            big_step = 0
        end if    

        if (flag=="c" .or. flag=="C") then ! checkpoint =====
            big_step = big_step + 1
            write (s,'(I4.4)') big_step
            call psi%save("sol_se" // s // ".h5")
            call hdf5_write_attributes("sol_se" // s // ".h5" , dnames = (/ "t" /), dvalues = (/ t /) )
        else
            call solution_out_schroedinger_1D(psi, t, dt, flag, abort_flag)
        end if   

        abort_flag = .false.
    end subroutine solution_out_se



    subroutine solution_out_gpe(psi, t, dt, flag, abort_flag)
        use tssm_hdf5
        class(wave_function), intent(inout) :: psi
        real(kind=prec), intent(in) :: t 
        real(kind=prec), intent(in) :: dt 
        character, intent(in) :: flag
        logical, intent(out) :: abort_flag
        integer, save :: big_step 
        character(len=4) :: s

        if (flag=="i" .or. flag=="I") then ! init =====
            big_step = 0
        end if    

        if (flag=="c" .or. flag=="C") then ! checkpoint =====
            big_step = big_step + 1
            write (s,'(I4.4)') big_step
            call psi%save("sol_gpe" // s // ".h5")
            call hdf5_write_attributes("sol_gpe" // s // ".h5" , dnames = (/ "t" /), dvalues = (/ t /) )
        else
            call solution_out_schroedinger_1D(psi, t, dt, flag, abort_flag)
        end if   

        abort_flag = .false.
    end subroutine solution_out_gpe

    subroutine solution_out_gpe_free(psi, t, dt, flag, abort_flag)
        use tssm_hdf5
        class(wave_function), intent(inout) :: psi
        real(kind=prec), intent(in) :: t 
        real(kind=prec), intent(in) :: dt 
        character, intent(in) :: flag
        logical, intent(out) :: abort_flag
        integer, save :: big_step 
        character(len=4) :: s

        if (flag=="i" .or. flag=="I") then ! init =====
            big_step = 0
        end if    

        if (flag=="c" .or. flag=="C") then ! checkpoint =====
            big_step = big_step + 1
            write (s,'(I4.4)') big_step
            call psi%save("sol_free" // s // ".h5")
            call hdf5_write_attributes("sol_free" // s // ".h5" , dnames = (/ "t" /), dvalues = (/ t /) )
        else
            call solution_out_schroedinger_1D(psi, t, dt, flag, abort_flag)
        end if   

        abort_flag = .false.
    end subroutine solution_out_gpe_free




  function harmonic_trap_1D(x) result(V)
        real(kind=prec), intent(in) :: x
        real(kind=prec) :: V
        V = 0.5_prec * x**2
  end function harmonic_trap_1D





end program simulation_1D
