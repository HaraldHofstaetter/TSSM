module gaussian_module
    use tssm_base, only: prec
    implicit none 

    real(prec) :: hbar = 1_prec
    real(prec) :: mass = 1_prec
    real(prec) :: a  = 3.0_prec
    real(prec) :: kx = +1.0_prec
    real(prec) :: ky = -0.5_prec
    real(prec) :: kz = 1.1_prec

    real(prec) :: t = 0_prec
    real(kind=prec), parameter :: pi = 4.0_prec*atan(1.0_prec)
contains
    function gaussian_1D(x)
        complex(kind=prec) :: gaussian_1D
        real(kind=prec), intent(in) :: x

        real(prec) :: phi, theta, a1, cx
        complex(prec) :: b
        !cf. C. Cohen-Tannoudji, B. Diu, F. Laloe: Quantenmechanik 2.Aufl., eq. (1.165)
        theta = 0.5_prec*atan(2.0_prec*hbar*t/(mass*a**2))
        phi = -theta - 0.5_prec*hbar*(kx**2)*t/mass
        a1 = (2.0_prec*a**2/pi/(a**4+4.0_prec*(hbar*t/mass)**2))**0.25_prec
        b = cmplx(a**2, 2.0_prec*hbar*t/mass, kind=prec)
        cx = hbar*kx*t/mass
        gaussian_1D = a1*exp(cmplx(0_prec, phi + kx*x, kind=prec) - ((x-cx)**2)/b)  
    end function gaussian_1D

    function gaussian_2D(x, y)
        complex(kind=prec) :: gaussian_2D
        real(kind=prec), intent(in) :: x, y

        real(prec) :: phi, theta, a1, cx, cy
        complex(prec) :: b


        !gaussian_2D = gaussian_1D(x)!*gaussian_1D(y)
        !return 

        theta = atan(2.0_prec*hbar*t/(mass*a**2))
        phi = -theta - 0.5_prec*hbar*(kx**2 + ky**2)*t/mass
        a1 = (2.0_prec*a**2/pi/(a**4+4.0_prec*(hbar*t/mass)**2))**0.5_prec
        b = cmplx(a**2, 2.0_prec*hbar*t/mass, kind=prec)
        cx = hbar*kx*t/mass
        cy = hbar*ky*t/mass
        gaussian_2D = a1*exp(cmplx(0_prec, phi + kx*x + ky*y, kind=prec) - ((x-cx)**2 + (y-cy)**2)/b)  
    end function gaussian_2D

    function gaussian_3D(x, y, z)
        complex(kind=prec) :: gaussian_3D
        real(kind=prec), intent(in) :: x, y, z

        real(prec) :: phi, theta, a1, cx, cy, cz
        complex(prec) :: b
        
        !gaussian_3D = gaussian_1D(x)*gaussian_1D(y)*gaussian_1D(z)
        !return 

        theta = 1.5_prec*atan(2.0_prec*hbar*t/(mass*a**2))
        phi = -theta - 0.5_prec*hbar*(kx**2 + ky**2 + kz**2)*t/mass
        a1 = (2.0_prec*a**2/pi/(a**4+4.0_prec*(hbar*t/mass)**2))**0.75_prec
        b = cmplx(a**2, 2.0_prec*hbar*t/mass, kind=prec)
        cx = hbar*kx*t/mass
        cy = hbar*ky*t/mass
        cz = hbar*kz*t/mass
       gaussian_3D = a1*exp(cmplx(0_prec, phi + kx*x + ky*y +kz*z, kind=prec) - ((x-cx)**2 + (y-cy)**2 + (z-cz)**2)/b)  
    end function gaussian_3D


end module gaussian_module


program test_gaussian
    use tssm_schroedinger
    use gaussian_module  
    implicit none 

    real(prec) :: tend, E_kin, E_pot, E, E_var 

    type(schroedinger_1D) :: m1, m11
    type(wf_schroedinger_1D) :: psi1, psi11
    type(wf_schroedinger_1D) :: psi1_ex

    type(schroedinger_2D) :: m2 
    type(wf_schroedinger_2D) :: psi2
    type(wf_schroedinger_2D) :: psi2_ex

    type(schroedinger_3D) :: m3 
    type(wf_schroedinger_3D) :: psi3
    type(wf_schroedinger_3D) :: psi3_ex


    call initialize_tssm

!------------------------------------------------------------
! 1D
!------------------------------------------------------------
#if 1
   m1 = schroedinger_1D(nx = 1024,           &
                       xmin = -32.0_prec,  &
                       xmax =  32.0_prec,  &
                       hbar = hbar,    &
                       mass = mass,    &
                       boundary_conditions = periodic)
    psi1 =  wf_schroedinger_1D(m1)
    psi1_ex =  wf_schroedinger_1D(m1)


    t = 0_prec
    call psi1%set(gaussian_1D) 
    call psi1%save("gaussian.h5") 
    call psi1%load("gaussian.h5") 

 !   call psi1%save("xxx.h5")
 !   m11 = schroedinger_1D(nx = 256,           &
 !                      xmin = -64.0_prec,  &
 !                      xmax =  64.0_prec,  &
 !                      hbar = hbar,    &
 !                      mass = mass,    &
 !                      boundary_conditions = periodic)
 !   psi11 =  wf_schroedinger_1D(m11)
 !
 !   call psi11%load("xxx.h5")
 !   call psi11%save("xxx1.h5")

    print *, "1D D", m1%g%dx
    print *, "1D SF", m1%g%n1min, m1%nf1min
    print *, "1D S1", shape(m1%eigenvalues1)
    print *, "1D S2", shape(psi1%uf)

    tend = 0.1_prec

    !get initial solution:
    t = 0_prec
    call psi1%set(gaussian_1D) 
    print *, "1D E_kin", psi1%kinetic_energy()
    call psi1%get_energy_expectation_deviation(E, E_var)
    print *, "1D E_kin, E_var", E, E_var

    !compute numerical solution
    call psi1%propagate_A(cmplx(tend,0,kind=prec))

    !get exact final solution
    t = tend 
    call psi1_ex%set(gaussian_1D) 

    print *, "1D E_kin", psi1%kinetic_energy()
    print *, "1D err", psi1%distance(psi1_ex)
#endif

!------------------------------------------------------------
! 2D
!------------------------------------------------------------
#if 0
   m2 = schroedinger_2D(nx = 128,           &
                  xmin = -30.0_prec,  &
                  xmax =  31.0_prec,  &
                  ny = 256,           &
                  ymin = -32.0_prec,  &
                  ymax =  33.0_prec,   &
                  hbar = hbar,    &
                  mass = mass,    &
                  boundary_conditions = dirichlet )

    psi2 =  wf_schroedinger_2D(m2)
    psi2_ex =  wf_schroedinger_2D(m2)

    tend = 1.0_prec

    !get initial solution:
    t = 0_prec
    call psi2%set(gaussian_2D) 
    print *, "2D  0 norm", psi2%norm()
    print *, "2D  0 E_kin", psi2%kinetic_energy()
    call psi2%get_energy_expectation_deviation(E, E_var)
    print *, "2D E_kin, E_var", E, E_var

    !compute numerical solution
    call psi2%propagate_A(cmplx(tend,0,kind=prec))
    print *, "2D P norm", psi2%norm()
    print *, "2D P E_kin", psi2%kinetic_energy()

    !get exact final solution
    t = tend 
    call psi2_ex%set(gaussian_2D) 
    print *, "2D E norm", psi2_ex%norm()
    print *, "2D E E_kin", psi2_ex%kinetic_energy()

    print *, "2D err", psi2%distance(psi2_ex)

    call psi2%save("xxx.h5")
    call psi2_ex%load("xxx.h5")
    call psi2_ex%save("xxx2.h5")

#endif    

#if 0
!------------------------------------------------------------
! 3D
!------------------------------------------------------------
  m3 = schroedinger_3D(nx = 128,           &
                  xmin = -30.0_prec,  &
                  xmax =  31.0_prec,  &
                  ny = 256,           &
                  ymin = -32.0_prec,  &
                  ymax =  33.0_prec,   &
                  nz = 128,           &
                  zmin = -34.0_prec,  &
                  zmax =  35.0_prec,   &
                  hbar = hbar,    &
                  mass = mass,    &
                  boundary_conditions = neumann)
    psi3 =  wf_schroedinger_3D(m3)
    psi3_ex =  wf_schroedinger_3D(m3)

    tend = 1.0_prec

    !get initial solution:
    t = 0_prec
    call psi3%set(gaussian_3D) 

    print *, "3D norm", psi3%norm()
    print *, "3D E_kin", psi3%kinetic_energy()
    call psi3%get_energy_expectation_deviation(E, E_var)
    print *, "3D E_kin, E_var", E, E_var
    
    !compute numerical solution
    call psi3%propagate_A(cmplx(tend,0,kind=prec))

    !get exact final solution
    t = tend 
    call psi3_ex%set(gaussian_3D) 

    print *, "3D norm", psi3%norm()
    print *, "3D err", psi3%distance(psi3_ex)

    call psi3%save("xxx.h5")
    call psi3_ex%load("xxx.h5")
    call psi3_ex%save("xxx2.h5")

#endif

    call finalize_tssm
end program test_gaussian
