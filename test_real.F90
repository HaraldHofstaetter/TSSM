

program test_real
   use tssm_fourier
   implicit none

   integer :: bc
   call initialize_tssm

   print *, "BC=PERIODIC"
   bc = periodic
#ifndef _MPI_
   call test_real_1D
#endif    
   call test_real_2D
   call test_real_3D
   print *, "BC=DIRICHLET"
   bc = dirichlet
#ifndef _MPI_   
   call test_real_1D
#endif    
   call test_real_2D
   call test_real_3D
   print *, "BC=NEUMANN"
   bc = neumann
#ifndef _MPI_   
   call test_real_1D
#endif    
   call test_real_2D
   call test_real_3D

   call finalize_tssm
contains   

subroutine test_real_1D
   type(fourier_1D) :: m
   type(wf_fourier_1D) :: psi

   type(fourier_real_1D) :: mr
   type(wf_fourier_real_1D) :: psir

   m = fourier_1D(nx=4, xmin=-2.0_prec, xmax=+1.0_prec, &
                  boundary_conditions = bc)
   mr = fourier_real_1D(nx=4, xmin=-2.0_prec, xmax=+1.0_prec, &
                        boundary_conditions = bc)
   psi = wf_fourier_1D(m)
   psir = wf_fourier_real_1D(mr)

   call psi%rset(f1)
   call psir%set(f1)

   call psi%propagate_A((.1_prec,0.0_prec))
   !call psi%to_frequency_space
   call psi%to_real_space

   call psir%propagate_A(.1_prec)
   !call psir%to_frequency_space
   call psir%to_real_space

   print *, "1D ERR", norm2(real(psi%u, kind=prec)-psir%u)
   !call psi%save("xxx0.h5") 
   !call psir%save("xxxr0.h5") 
end subroutine test_real_1D

subroutine test_real_2D
   type(fourier_2D) :: m
   type(wf_fourier_2D) :: psi

   type(fourier_real_2D) :: mr
   type(wf_fourier_real_2D) :: psir

   m = fourier_2D(nx=6, xmin=-2.0_prec, xmax=+1.0_prec, &
                  ny=18, ymin=-3.0_prec, ymax=+1.0_prec, &
                  boundary_conditions = bc)
   mr = fourier_real_2D(nx=6, xmin=-2.0_prec, xmax=+1.0_prec, &
                        ny=18, ymin=-3.0_prec, ymax=+1.0_prec, &
                        boundary_conditions = bc)
   psi = wf_fourier_2D(m)
   psir = wf_fourier_real_2D(mr)

   call psi%rset(f2)
   call psir%set(f2)

   call psi%propagate_A((.1_prec,0.0_prec))
   call psi%to_real_space

   call psir%propagate_A(.1_prec)
   call psir%to_real_space

   print *, "2D ERR", norm2(real(psi%u, kind=prec)-psir%u)
end subroutine test_real_2D

subroutine test_real_3D
   type(fourier_3D) :: m
   type(wf_fourier_3D) :: psi

   type(fourier_real_3D) :: mr
   type(wf_fourier_real_3D) :: psir

   m = fourier_3D(nx=4, xmin=-2.0_prec, xmax=+1.0_prec, &
                  ny=8, ymin=-3.0_prec, ymax=+1.0_prec, &
                  nz=12, zmin=-1.0_prec, zmax=+1.0_prec, &
                  boundary_conditions = bc)
   mr = fourier_real_3D(nx=4, xmin=-2.0_prec, xmax=+1.0_prec, &
                        ny=8, ymin=-3.0_prec, ymax=+1.0_prec, &
                        nz=12, zmin=-1.0_prec, zmax=+1.0_prec, &
                        boundary_conditions = bc)
   psi = wf_fourier_3D(m)
   psir = wf_fourier_real_3D(mr)

   call psi%rset(f3)
   call psir%set(f3)

   call psi%propagate_A((.1_prec,0.0_prec))
!   call psi%to_frequency_space
   call psi%to_real_space

   call psir%propagate_A(.1_prec)
!   call psir%to_frequency_space
   call psir%to_real_space

   print *, "3D ERR", norm2(real(psi%u, kind=prec)-psir%u)
end subroutine test_real_3D

   function f1(x)
   real(prec):: f1
   real(prec), intent(in) :: x
   f1 = x**2
   end function f1

   function f2(x,y)
   real(prec):: f2
   real(prec), intent(in) :: x,y
   f2 = x**2+y**2
   end function f2

   function f3(x,y,z)
   real(prec):: f3
   real(prec), intent(in) :: x,y,z
   f3 = x**2+y**2+z**2
   end function f3

   subroutine real_mat_out(A)
    real(kind=8), intent(in) :: A(:,:)

    integer :: i,j
    do i=1,ubound(A,1)
        !write(*,91)(A(i,j),j=1,n)
        write(*,91)(real(A(i,j), kind=8),j=1,ubound(A,2))
    end do
 91 format(1x,2049(E11.5,1x))    
end subroutine real_mat_out

end program test_real
