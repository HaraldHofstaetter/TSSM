program test_imag_time   
    use tssm_schroedinger
    implicit none

   integer :: bc
   call initialize_tssm

   if (this_proc==0) then 
      print *, "BC=PERIODIC"
   end if    
   bc = periodic
#ifndef _MPI_
   call test_real_1D
#endif    
   call test_real_2D
   call test_real_3D
   if (this_proc==0) then 
       print *, "BC=DIRICHLET"
   end if    
   bc = dirichlet
#ifndef _MPI_   
   call test_real_1D
#endif    
   call test_real_2D
   call test_real_3D
   if (this_proc==0) then 
       print *, "BC=NEUMANN"
   end if    
   bc = neumann
#ifndef _MPI_   
   call test_real_1D
#endif    
   call test_real_2D
   call test_real_3D

   call finalize_tssm
contains   

subroutine test_real_1D
   real(prec) :: E_exp, E_dev, E_kin, E_int, E_pot, obs, N, Nf, x_mean, x_dev
   type(schroedinger_1D) :: m
   type(wf_schroedinger_1D) :: psi, psi1

   type(schroedinger_real_1D) :: mr
   type(wf_schroedinger_real_1D) :: psir, psi1r


   m = schroedinger_1D(nx=32, xmin=-1.0_prec, xmax=+1.0_prec, &
                  potential = h1, &
                  boundary_conditions = bc)
   mr = schroedinger_real_1D(nx=32, xmin=-1.0_prec, xmax=+1.0_prec, &
                        potential = h1, &
                        boundary_conditions = bc)
   psi = wf_schroedinger_1D(m)
   psir = wf_schroedinger_real_1D(mr)

   psi1 = wf_schroedinger_1D(m)
   psi1r = wf_schroedinger_real_1D(mr)

   call psi%rset(f1)
   psi1%u = 0.0_prec
   call psi%add_apply_A(psi1)
   call psi%get_energy_expectation_deviation(E_exp, E_dev)
   E_kin = psi%kinetic_energy()
   E_pot = psi%potential_energy()
   obs = psi%observable(h1)
   N = psi%norm()
   Nf = psi%norm_in_frequency_space()
   if (this_proc==0) then 
       print *, "1D E_KIN full", E_kin+E_pot, E_exp, E_dev, E_pot, obs, N, Nf
   end if    
   call psi%get_realspace_observables(E_pot, E_int, x_mean, x_dev)
   if (this_proc==0) then 
       print *, "1D OBS full", E_pot, E_int, x_mean, x_dev
   end if    

   call psir%set(f1)
   psi1r%u = 0.0_prec
   call psir%add_apply_A(psi1r)
   call psir%get_energy_expectation_deviation(E_exp, E_dev)
   E_kin = psir%kinetic_energy()
   E_pot = psir%potential_energy()
   obs = psir%observable(h1)
   N = psir%norm()
   Nf = psir%norm_in_frequency_space()
   if (this_proc==0) then 
       print *, "1D E_KIN real", E_kin+E_pot, E_exp, E_dev, E_pot, obs, N, Nf
   end if    
   call psir%get_realspace_observables(E_pot, E_int, x_mean, x_dev)
   if (this_proc==0) then 
       print *, "1D OBS real", E_pot, E_int, x_mean, x_dev**2
   end if    
end subroutine test_real_1D

function f1(x)
   real(prec):: f1
   real(prec), intent(in) :: x
   f1 = cos(2*pi*x)
   !f1 = sin(2*pi*x)
end function f1

function h1(x)
   real(prec):: h1
   real(prec), intent(in) :: x
   h1 = x**2
end function h1




subroutine test_real_2D
   real(prec) :: E_exp, E_dev, E_int, E_kin, E_pot, obs, N, Nf, x_mean, x_dev, y_mean, y_dev
   type(schroedinger_2D) :: m
   type(wf_schroedinger_2D) :: psi, psi1

   type(schroedinger_real_2D) :: mr
   type(wf_schroedinger_real_2D) :: psir, psi1r

   m = schroedinger_2D(nx=32, xmin=-1.0_prec, xmax=+1.0_prec, &
                       ny=32, ymin=-1.0_prec, ymax=+1.0_prec, &
                       potential = h2, &
                       boundary_conditions = bc)
   mr = schroedinger_real_2D(nx=32, xmin=-1.0_prec, xmax=+1.0_prec, &
                                       ny=32, ymin=-1.0_prec, ymax=+1.0_prec, &
                                       potential = h2, &
                                       boundary_conditions = bc)
   psi = wf_schroedinger_2D(m)
   psir = wf_schroedinger_real_2D(mr)

   psi1 = wf_schroedinger_2D(m)
   psi1r = wf_schroedinger_real_2D(mr)

   call psi%rset(f2)
   psi1%u = 0.0_prec
   call psi%add_apply_A(psi1)
   call psi%get_energy_expectation_deviation(E_exp, E_dev)
   E_kin = psi%kinetic_energy()
   E_pot = psi%potential_energy()
   obs = psi%observable(h2)
   N = psi%norm()
   Nf = psi%norm_in_frequency_space()
   if (this_proc==0) then 
       print *, "2D E_KIN full", E_kin+E_pot, E_exp, E_dev, E_pot, obs, N, Nf
   end if    
   call psi%get_realspace_observables(E_pot, E_int, x_mean, x_dev,y_mean, y_dev)
   if (this_proc==0) then 
       print *, "2D OBS full", E_pot, E_int, x_mean, y_mean, x_dev**2+y_dev**2
   end if    

   call psir%set(f2)
   psi1r%u = 0.0_prec
   call psir%add_apply_A(psi1r)
   call psir%get_energy_expectation_deviation(E_exp, E_dev)
   E_kin = psir%kinetic_energy()
   E_pot = psir%potential_energy()
   obs = psir%observable(h2)
   N = psir%norm()
   Nf = psir%norm_in_frequency_space()
   if (this_proc==0) then 
       print *, "2D E_KIN real", E_kin+E_pot, E_exp, E_dev, E_pot, obs, N, Nf
   end if    
   call psir%get_realspace_observables(E_pot, E_int, x_mean, x_dev,y_mean, y_dev)
   if (this_proc==0) then 
       print *, "2D OBS full", E_pot, E_int, x_mean, y_mean, x_dev**2+y_dev**2
   end if    
end subroutine test_real_2D




function f2(x, y)
   real(prec):: f2
   real(prec), intent(in) :: x,y
   !f2 = cos(2*pi*x)*cos(2*pi*y)
   f2 = sin(2*pi*x)*sin(2*pi*y)
   !f2 = x**2-3*y**2
end function f2

function h2(x,y)
   real(prec):: h2
   real(prec), intent(in) :: x, y
   h2 = x**2 + y**2
end function h2




subroutine test_real_3D
   real(prec) :: E_exp, E_dev, E_kin, E_pot, E_int, obs, N, Nf, x_mean, x_dev,y_mean, y_dev,z_mean, z_dev
   type(schroedinger_3D) :: m
   type(wf_schroedinger_3D) :: psi, psi1

   type(schroedinger_real_3D) :: mr
   type(wf_schroedinger_real_3D) :: psir, psi1r

   m = schroedinger_3D(nx=32, xmin=-1.0_prec, xmax=+1.0_prec, &
                       ny=32, ymin=-1.0_prec, ymax=+1.0_prec, &
                       nz=32, zmin=-1.0_prec, zmax=+1.0_prec, &
                       potential = h3, &
                       boundary_conditions = bc)
   mr = schroedinger_real_3D(nx=32, xmin=-1.0_prec, xmax=+1.0_prec, &
                                       ny=32, ymin=-1.0_prec, ymax=+1.0_prec, &
                                       nz=32, zmin=-1.0_prec, zmax=+1.0_prec, &
                                       potential = h3, &
                                       boundary_conditions = bc)
   psi = wf_schroedinger_3D(m)
   psir = wf_schroedinger_real_3D(mr)

   psi1 = wf_schroedinger_3D(m)
   psi1r = wf_schroedinger_real_3D(mr)

   call psi%rset(f3)
   psi1%u = 0.0_prec
   call psi%add_apply_A(psi1)
   call psi%get_energy_expectation_deviation(E_exp, E_dev)
   E_kin = psi%kinetic_energy()
   E_pot = psi%potential_energy()
   obs = psi%observable(h3)
   N = psi%norm()
   Nf = psi%norm_in_frequency_space()
   if (this_proc==0) then 
        print *, "3D E_KIN full", E_kin+E_pot, E_exp, E_dev, E_pot, obs, N, Nf
   end if    
   call psi%get_realspace_observables(E_pot, E_int, x_mean, x_dev,y_mean, y_dev, z_mean, z_dev)
   if (this_proc==0) then 
       print *, "3D OBS full", E_pot, E_int, x_mean, y_mean,  z_mean,x_dev**2+y_dev**2+z_dev**2
   end if    

   call psir%set(f3)
   psi1r%u = 0.0_prec
   call psir%add_apply_A(psi1r)
   call psir%get_energy_expectation_deviation(E_exp, E_dev)
   E_kin = psir%kinetic_energy()
   E_pot = psir%potential_energy()
   obs = psir%observable(h3)
   N = psir%norm()
   Nf = psir%norm_in_frequency_space()
   if (this_proc==0) then 
        print *, "3D E_KIN real", E_kin+E_pot, E_exp, E_dev, E_pot, obs, N, Nf
   end if    
   call psir%get_realspace_observables(E_pot, E_int, x_mean, x_dev,y_mean, y_dev, z_mean, z_dev)
   if (this_proc==0) then 
       print *, "3D OBS real", E_pot, E_int, x_mean, y_mean, z_mean, x_dev**2+y_dev**2+z_dev**2
   end if    
end subroutine test_real_3D




function f3(x, y, z)
   real(prec):: f3
   real(prec), intent(in) :: x,y, z
   !f3 = cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z)
   f3 = sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)
   !f3 = x**2-3*y**2+2*z**2
end function f3

function h3(x,y,z)
   real(prec):: h3
   real(prec), intent(in) :: x, y, z
   h3 = x**2 + y**2 + z**2
end function h3

   


end program test_imag_time
