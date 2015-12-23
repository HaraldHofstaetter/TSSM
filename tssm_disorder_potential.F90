module tssm_disorder_potential
    use tssm_schroedinger
    implicit none 


contains

    subroutine compute_disorder_potential_1D(g, V, dx, sigma, f) 
        class(grid_equidistant_1D) :: g
        real(kind=prec), intent(inout) :: V(g%n1min:g%n1max)
        real(kind=prec), intent(in) :: dx, sigma
        integer, intent(in) :: f

        integer :: shell, shells
        integer :: ix
        real(prec) :: a, cx

        V = 0.0_prec

        shells = floor(maxval( (/ abs(g%xmin-f*sigma)/dx, abs(g%xmax+f*sigma)/dx /) ))

        call random_number(a)
        call add_gaussian_1D(g, V, a-0.5_prec, 0.0_prec, sigma, f)

        do shell = 1, shells
            do ix = -shell,+shell,2*shell
                cx = ix*dx
                call random_number(a)
                call add_gaussian_1D(g, V, a-0.5_prec, cx, sigma, f)
            end do
        end do

    end subroutine compute_disorder_potential_1D


    subroutine compute_disorder_potential_2D(g, V, dx, dy, sigma, f) 
        class(grid_equidistant_2D) :: g
        real(kind=prec), intent(inout) :: V(g%n1min:g%n1max, g%n2min:g%n2max)
        real(kind=prec), intent(in) :: dx, dy, sigma
        integer, intent(in) :: f

        integer :: shell, shells
        integer :: ix, iy, k
        real(prec) :: a, cx, cy

        V = 0.0_prec

        shells = floor(maxval( (/ abs(g%xmin-f*sigma)/dx, abs(g%xmax+f*sigma)/dx, &
                                  abs(g%ymin-f*sigma)/dy, abs(g%ymax+f*sigma)/dy /) ))

        call random_number(a)
        call add_gaussian_2D(g, V, a-0.5_prec, 0.0_prec, 0.0_prec, sigma, f)
        
        do shell = 1, shells
           k=0
           do iy = -shell,shell, 2*shell
                cy = iy*dy
                do ix = -shell, shell
                    cx = ix*dx
                    call random_number(a)
                    call add_gaussian_2D(g, V, a-0.5_prec, cx, cy, sigma, f)
                    k = k+1
                end do
            end do

            do ix = -shell,shell, 2*shell
                cx = ix*dx
                do iy = -shell+1, shell-1
                    cy = iy*dy
                    call random_number(a)
                    call add_gaussian_2D(g, V, a-0.5_prec, cx, cy, sigma, f)
                    k = k+1
                end do
            end do
        
           print *, "shell", shell, shells, k
        end do



    end subroutine compute_disorder_potential_2D



    subroutine compute_disorder_potential_3D(g, V, dx, dy, dz, sigma, f) 
        class(grid_equidistant_3D) :: g
        real(kind=prec), intent(inout) :: V(g%n1min:g%n1max, g%n2min:g%n2max, g%n3min:g%n3max)
        real(kind=prec), intent(in) :: dx, dy, dz, sigma
        integer, intent(in) :: f

        integer :: shell, shells
        integer :: ix, iy, iz
        real(prec) :: a, cx, cy, cz

        shells = floor(maxval( (/ abs(g%xmin-f*sigma)/dx, abs(g%xmax+f*sigma)/dx, &
                                  abs(g%ymin-f*sigma)/dy, abs(g%ymax+f*sigma)/dy, &
                                  abs(g%zmin-f*sigma)/dz, abs(g%zmax+f*sigma)/dz /) ))

        V = 0.0_prec

        call random_number(a)
        call add_gaussian_3D(g, V, a-0.5_prec, 0.0_prec, 0.0_prec, 0.0_prec, sigma, f)

        do shell = 1, shells
            do iz = -shell,shell, 2*shell
                cz = iz*dz
                do iy = -shell,shell
                    cy = iy*dy
                    do ix = -shell, shell
                        cx = ix*dx
                        call random_number(a)
                        call add_gaussian_3D(g, V, a-0.5_prec, cx, cy, cz, sigma, f)
                    end do
                end do
            end do

            do iy = -shell,shell, 2*shell
                cy = iy*dy
                do iz = -shell+1,shell-1
                    cz = iz*dz
                    do ix = -shell, shell
                        cx = ix*dx
                        call random_number(a)
                        call add_gaussian_3D(g, V, a-0.5_prec, cx, cy, cz, sigma, f)
                    end do
                end do
            end do

            do ix = -shell,shell, 2*shell
                cx = ix*dx
                do iz = -shell+1,shell-1
                    cz = iz*dz
                    do iy = -shell+1, shell-1
                        cy = iy*dy
                        call random_number(a)
                        call add_gaussian_3D(g, V, a-0.5_prec, cx, cy, cz, sigma, f)
                    end do
                end do
            end do

        end do
        


    end subroutine compute_disorder_potential_3D





    subroutine add_gaussian_1D(g, V, a, cx, sigma, f)
        class(grid_equidistant_1D) :: g
        real(kind=prec), intent(inout) :: V(g%n1min:g%n1max)
        real(kind=prec), intent(in) :: a, cx, sigma
        integer, intent(in) :: f

        integer :: n1min, n1max
        integer :: ix
        real(prec) :: x

         n1min = max(g%n1min, floor((cx-g%xmin-f*sigma)/g%dx))
         n1max = min(g%n1max, ceiling((cx-g%xmin+f*sigma)/g%dx))


         do ix=n1min, n1max
             x = g%xmin + ix*g%dx
             V(ix) = V(ix) + a*exp(-((x-cx)**2)/(2.0_prec*sigma**2))
         end do    

    end subroutine add_gaussian_1D


    subroutine add_gaussian_2D(g, V, a, cx, cy, sigma, f)
        class(grid_equidistant_2D) :: g
        real(kind=prec), intent(inout) :: V(g%n1min:g%n1max, g%n2min:g%n2max)
        real(kind=prec), intent(in) :: a, cx, cy, sigma
        integer, intent(in) :: f

        integer :: n1min, n1max, n2min, n2max
        integer :: ix, iy
        real(prec) :: x, y 

         n1min = max(g%n1min, floor((cx-g%xmin-f*sigma)/g%dx))
         n1max = min(g%n1max, ceiling((cx-g%xmin+f*sigma)/g%dx))
         n2min = max(g%n2min, floor((cy-g%ymin-f*sigma)/g%dy))
         n2max = min(g%n2max, ceiling((cy-g%ymin+f*sigma)/g%dy))


        do iy = n2min, n2max
            y = g%ymin + iy*g%dy
            do ix=n1min, n1max
                x = g%xmin + ix*g%dx
                V(ix,iy) = V(ix,iy) + a*exp(-((x-cx)**2+(y-cy)**2)/(2.0_prec*sigma**2))
            end do    
        end do    

    end subroutine add_gaussian_2D


   subroutine add_gaussian_3D(g, V, a, cx, cy, cz, sigma, f)
        class(grid_equidistant_3D) :: g
        real(kind=prec), intent(inout) :: V(g%n1min:g%n1max, g%n2min:g%n2max, g%n3min:g%n3max)
        real(kind=prec), intent(in) :: a, cx, cy, cz, sigma
        integer, intent(in) :: f

        integer :: n1min, n1max, n2min, n2max, n3min, n3max
        integer :: ix, iy, iz
        real(prec) :: x, y, z 

         n1min = max(g%n1min, floor((cx-g%xmin-f*sigma)/g%dx))
         n1max = min(g%n1max, ceiling((cx-g%xmin+f*sigma)/g%dx))
         n2min = max(g%n2min, floor((cy-g%ymin-f*sigma)/g%dy))
         n2max = min(g%n2max, ceiling((cy-g%ymin+f*sigma)/g%dy))
         n3min = max(g%n3min, floor((cz-g%zmin-f*sigma)/g%dz))
         n3max = min(g%n3max, ceiling((cz-g%zmin+f*sigma)/g%dz))


        do iz = n3min, n3max
            z = g%zmin + iz*g%dz
            do iy = n2min, n2max
                y = g%ymin + iy*g%dy
                do ix=n1min, n1max
                    x = g%xmin + ix*g%dx
                    V(ix,iy,iz) = V(ix,iy,iz) + a*exp(-((x-cx)**2+(y-cy)**2+(z-cz)**2)&
                                                                 /(2.0_prec*sigma**2))
                end do    
            end do    
        end do    

    end subroutine add_gaussian_3D



    subroutine get_mean_value_and_standard_deviation_1D(g, V, mean_value, standard_deviation)
        class(grid_equidistant_1D) :: g
        real(kind=prec), intent(inout) :: V(g%n1min:g%n1max)   
        real(kind=prec), intent(out) :: mean_value, standard_deviation
        
        real(kind=prec) :: s, s1
#ifdef _MPI_     
        integer :: ierr
#endif        

        s = sum(V)

#ifdef _MPI_     
        call MPI_Reduce(s, s1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            mean_value = s1/(g%nn1max-g%nn1min+1)
        end if    
        call MPI_Bcast(mean_value, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#else
         mean_value = s/(g%nn1max-g%nn1min+1)
#endif 

         s = sum((V-mean_value)**2)

#ifdef _MPI_     
        call MPI_Reduce(s, s1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            standard_deviation = sqrt(s1/(g%nn1max-g%nn1min+1))
        end if    
        call MPI_Bcast(standard_deviation, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#else
         standard_deviation = sqrt(s/(g%nn1max-g%nn1min+1))
#endif 
    end subroutine get_mean_value_and_standard_deviation_1D
    


    subroutine get_mean_value_and_standard_deviation_2D(g, V, mean_value, standard_deviation)
        class(grid_equidistant_2D) :: g
        real(kind=prec), intent(inout) :: V(g%n1min:g%n1max, g%n2min:g%n2max)   
        real(kind=prec), intent(out) :: mean_value, standard_deviation
        
        real(kind=prec) :: s, s1
#ifdef _MPI_     
        integer :: ierr
#endif        

        s = sum(V)

#ifdef _MPI_     
        call MPI_Reduce(s, s1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            mean_value = s1/((g%nn1max-g%nn1min+1)*(g%nn2max-g%nn2min+1))
        end if    
        call MPI_Bcast(mean_value, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#else
         mean_value = s/((g%nn1max-g%nn1min+1)*(g%nn2max-g%nn2min+1))
#endif 

         s = sum((V-mean_value)**2)

#ifdef _MPI_     
        call MPI_Reduce(s, s1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            standard_deviation = sqrt(s1/((g%nn1max-g%nn1min+1)*(g%nn2max-g%nn2min+1)))
        end if    
        call MPI_Bcast(standard_deviation, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#else
        standard_deviation = sqrt(s/((g%nn1max-g%nn1min+1)*(g%nn2max-g%nn2min+1)))
#endif 
    end subroutine get_mean_value_and_standard_deviation_2D



    subroutine get_mean_value_and_standard_deviation_3D(g, V, mean_value, standard_deviation)
        class(grid_equidistant_3D) :: g
        real(kind=prec), intent(inout) :: V(g%n1min:g%n1max, g%n2min:g%n2max, g%n3min:g%n3max)   
        real(kind=prec), intent(out) :: mean_value, standard_deviation
        
        real(kind=prec) :: s, s1
#ifdef _MPI_     
        integer :: ierr
#endif        

        s = sum(V)

#ifdef _MPI_     
        call MPI_Reduce(s, s1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            mean_value = s1/((g%nn1max-g%nn1min+1)*(g%nn2max-g%nn2min+1)*(g%nn3max-g%nn3min+1))
        end if    
        call MPI_Bcast(mean_value, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#else
         mean_value = s/((g%nn1max-g%nn1min+1)*(g%nn2max-g%nn2min+1)*(g%nn3max-g%nn3min+1))
#endif 

         s = sum((V-mean_value)**2)

#ifdef _MPI_     
        call MPI_Reduce(s, s1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            standard_deviation = sqrt(s1/((g%nn1max-g%nn1min+1)*(g%nn2max-g%nn2min+1)*(g%nn3max-g%nn3min+1)))
        end if    
        call MPI_Bcast(standard_deviation, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#else
         standard_deviation = sqrt(s/((g%nn1max-g%nn1min+1)*(g%nn2max-g%nn2min+1)*(g%nn3max-g%nn3min+1)))
#endif 
    end subroutine get_mean_value_and_standard_deviation_3D

           


end module tssm_disorder_potential
