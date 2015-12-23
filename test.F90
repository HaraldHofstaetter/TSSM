program test
    use tssm
    use tssm_hdf5
    use tssm_fourier_1D
    implicit none

    character(len=20) :: gridname
    integer :: len

    type(fourier_1D) :: m, m1 
    type(wf_fourier_1D) :: psi, psi1

    type(grid_equidistant_1D) :: g

    m = fourier_1D(nx=128, xmin=-4.0_prec, xmax=4.0_prec, boundary_conditions=periodic) !,  ny=8, ymin=-5.0d0, ymax=5.0d0)
    m1 = fourier_1D(nx=128, xmin=-4.0_prec, xmax=4.0_prec, boundary_conditions=dirichlet) !,  ny=8, ymin=-5.0d0, ymax=5.0d0)
    psi =  wf_fourier_1D(m)
    psi1 =  wf_fourier_1D(m1)

    call psi%rset(f)
    call psi%save('1.h5')

   call hdf5_read_string_attribute("1.h5", "grid", gridname, len) 
   print *, "len",len 
   print *, "GRIDNAME", gridname(1:len)
   call hdf5_read_grid_attributes(g, "1.h5", "psi_real") 

   print *, "nx", g%nx, m%g%nx
   print *, "xmin", g%xmin, m%g%xmin
   print *, "xmax", g%xmax, m%g%xmax  
   print *, "nn1min", g%nn1min, m%g%nn1min
   print *, "nn1max", g%nn1max,  m%g%nn1max

    call psi1%load('1.h5')
    call psi1%save('2.h5')
    


contains
   function f(x)
       real(prec) :: f
       real(prec), intent(in) :: x
       f=cos(pi/4_prec*x)
   end function f

   function df(x)
       real(prec) :: df
       real(prec), intent(in) :: x
       df=-(pi/4_prec)**2*cos(pi/4_prec*x)
   end function df


   function gaussian(x ,y)
       real(prec) :: gaussian
       real(prec), intent(in) :: x, y
       gaussian = exp(-(x**2+y**2))
   end function gaussian



end program test
