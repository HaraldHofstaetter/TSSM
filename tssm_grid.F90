#ifdef _QUADPRECISION_
module tssmq_grid
    use tssmq_common
#else
module tssm_grid
    use tssm_common
#endif    
    use, intrinsic :: iso_c_binding, only: c_size_t
    implicit none

    type, abstract :: grid
        integer(kind=C_SIZE_T) :: alloc_size
#ifdef _OPENMP
        integer, allocatable :: jj(:)
#endif
    contains
    end type grid

    type, abstract, extends(grid) :: grid_1D
        integer :: nn1min, nn1max ! overall index range
        integer :: n1min, n1max ! this index section "belongs" to current processor
        integer :: m1min, m1max ! this index section (superset of the above) of a grid function
                                ! is allocated by current processor
    contains
        procedure :: allocate_real_gridfun => allocate_real_gridfun_1D
        procedure :: allocate_complex_gridfun => allocate_complex_gridfun_1D
        procedure :: scale_real_gridfun => scale_real_gridfun_1D
        procedure :: scale_complex_gridfun => scale_complex_gridfun_1D
        procedure :: axpy_real_gridfun => axpy_real_gridfun_1D
        procedure :: axpy_complex_gridfun => axpy_complex_gridfun_1D
        procedure :: copy_real_gridfun => copy_real_gridfun_1D
        procedure :: copy_complex_gridfun => copy_complex_gridfun_1D
    end type grid_1D

    type, abstract, extends(grid) :: grid_2D
        integer :: nn1min, nn1max 
        integer :: n1min, n1max
        integer :: m1min, m1max
        integer :: nn2min, nn2max 
        integer :: n2min, n2max
        integer :: m2min, m2max
    contains
        procedure :: allocate_real_gridfun => allocate_real_gridfun_2D
        procedure :: allocate_complex_gridfun => allocate_complex_gridfun_2D
        procedure :: scale_real_gridfun => scale_real_gridfun_2D
        procedure :: axpy_real_gridfun => axpy_real_gridfun_2D
        procedure :: copy_real_gridfun => copy_real_gridfun_2D
        procedure :: scale_complex_gridfun => scale_complex_gridfun_2D
        procedure :: axpy_complex_gridfun => axpy_complex_gridfun_2D
        procedure :: copy_complex_gridfun => copy_complex_gridfun_2D
    end type grid_2D

    type, abstract, extends(grid) :: grid_3D
        integer :: nn1min, nn1max 
        integer :: n1min, n1max
        integer :: m1min, m1max
        integer :: nn2min, nn2max 
        integer :: n2min, n2max
        integer :: m2min, m2max
        integer :: nn3min, nn3max 
        integer :: n3min, n3max
        integer :: m3min, m3max
    contains
        procedure :: allocate_real_gridfun => allocate_real_gridfun_3D
        procedure :: allocate_complex_gridfun => allocate_complex_gridfun_3D
        procedure :: scale_real_gridfun => scale_real_gridfun_3D
        procedure :: axpy_real_gridfun => axpy_real_gridfun_3D
        procedure :: copy_real_gridfun => copy_real_gridfun_3D
        procedure :: scale_complex_gridfun => scale_complex_gridfun_3D
        procedure :: axpy_complex_gridfun => axpy_complex_gridfun_3D
        procedure :: copy_complex_gridfun => copy_complex_gridfun_3D
    end type grid_3D

!    interface allocate_gridfun_1D
!       module procedure  allocate_real_gridfun_1D
!       module procedure  allocate_complex_gridfun_1D
!    end interface allocate_gridfun_1D
!
!    interface allocate_gridfun_2D
!       module procedure  allocate_real_gridfun_2D
!       module procedure  allocate_complex_gridfun_2D
!    end interface allocate_gridfun_2D
!
!    interface allocate_gridfun_3D
!       module procedure  allocate_real_gridfun_3D
!       module procedure  allocate_complex_gridfun_3D
!    end interface allocate_gridfun_3D


    type, extends(grid_1D) :: grid_equidistant_1D
        integer :: nx
        real (kind=prec) :: xmin
        real (kind=prec) :: xmax
        real (kind=prec) :: dx
        real(kind=prec), allocatable :: nodes_x(:)
    contains
        procedure :: set_real_gridfun => set_real_gridfun_equidistant_1D
        procedure :: set_complex_gridfun => set_complex_gridfun_equidistant_1D
        procedure :: rset_complex_gridfun => rset_complex_gridfun_equidistant_1D
        procedure :: set_t_real_gridfun => set_t_real_gridfun_equidistant_1D
        procedure :: set_t_complex_gridfun => set_t_complex_gridfun_equidistant_1D
        procedure :: rset_t_complex_gridfun => rset_t_complex_gridfun_equidistant_1D
        procedure :: norm_real_gridfun => norm_real_gridfun_equidistant_1D
        procedure :: norm_complex_gridfun => norm_complex_gridfun_equidistant_1D
        procedure :: inner_product_real_gridfun => inner_product_real_gridfun_equidistant_1D
        procedure :: inner_product_complex_gridfun => inner_product_complex_gridfun_equidistant_1D
    end type grid_equidistant_1D

    interface grid_equidistant_1D
        module procedure new_grid_equidistant_1D 
    end interface grid_equidistant_1D

    type, extends(grid_2D) :: grid_equidistant_2D
        integer :: nx
        real (kind=prec) :: xmin
        real (kind=prec) :: xmax
        real (kind=prec) :: dx
        real(kind=prec), allocatable :: nodes_x(:)
        integer :: ny
        real (kind=prec) :: ymin
        real (kind=prec) :: ymax
        real (kind=prec) :: dy
        real(kind=prec), allocatable :: nodes_y(:)
    contains
        procedure :: set_real_gridfun => set_real_gridfun_equidistant_2D
        procedure :: set_complex_gridfun => set_complex_gridfun_equidistant_2D
        procedure :: rset_complex_gridfun => rset_complex_gridfun_equidistant_2D
        procedure :: set_t_real_gridfun => set_t_real_gridfun_equidistant_2D
        procedure :: set_t_complex_gridfun => set_t_complex_gridfun_equidistant_2D
        procedure :: rset_t_complex_gridfun => rset_t_complex_gridfun_equidistant_2D
        procedure :: norm_real_gridfun => norm_real_gridfun_equidistant_2D
        procedure :: norm_complex_gridfun => norm_complex_gridfun_equidistant_2D
!        procedure :: inner_product_real_gridfun => inner_product_real_gridfun_equidistant_2D
!        procedure :: inner_product_complex_gridfun => inner_product_complex_gridfun_equidistant_2D
    end type grid_equidistant_2D

    interface grid_equidistant_2D
        module procedure new_grid_equidistant_2D 
    end interface grid_equidistant_2D

    type, extends(grid_3D) :: grid_equidistant_3D
        integer :: nx
        real (kind=prec) :: xmin
        real (kind=prec) :: xmax
        real (kind=prec) :: dx
        real(kind=prec), allocatable :: nodes_x(:)
        integer :: ny
        real (kind=prec) :: ymin
        real (kind=prec) :: ymax
        real (kind=prec) :: dy
        real(kind=prec), allocatable :: nodes_y(:)
        integer :: nz
        real (kind=prec) :: zmin
        real (kind=prec) :: zmax
        real (kind=prec) :: dz
        real(kind=prec), allocatable :: nodes_z(:)
    contains
        procedure :: set_real_gridfun => set_real_gridfun_equidistant_3D
        procedure :: set_complex_gridfun => set_complex_gridfun_equidistant_3D
        procedure :: rset_complex_gridfun => rset_complex_gridfun_equidistant_3D
        procedure :: set_t_real_gridfun => set_t_real_gridfun_equidistant_3D
        procedure :: set_t_complex_gridfun => set_t_complex_gridfun_equidistant_3D
        procedure :: rset_t_complex_gridfun => rset_t_complex_gridfun_equidistant_3D
        procedure :: norm_real_gridfun => norm_real_gridfun_equidistant_3D
        procedure :: norm_complex_gridfun => norm_complex_gridfun_equidistant_3D
!        procedure :: inner_productinner_product_real_gridfun => inner_product_real_gridfun_equidistant_3D
!        procedure :: inner_product_complex_gridfun => inner_product_complex_gridfun_equidistant_3D
    end type grid_equidistant_3D

    interface grid_equidistant_3D
        module procedure new_grid_equidistant_3D 
    end interface grid_equidistant_3D


    type, extends(grid_1D) :: grid_tensorial_1D
        integer :: nx
        real(kind=prec), allocatable :: nodes_x(:)
        !real(kind=prec), allocatable :: weights_x(:)
        real(kind=prec), pointer :: weights_x(:)
    contains
        procedure :: set_real_gridfun => set_real_gridfun_tensorial_1D
        procedure :: set_complex_gridfun => set_complex_gridfun_tensorial_1D
        procedure :: rset_complex_gridfun => rset_complex_gridfun_tensorial_1D
        procedure :: set_t_real_gridfun => set_t_real_gridfun_tensorial_1D
        procedure :: set_t_complex_gridfun => set_t_complex_gridfun_tensorial_1D
        procedure :: rset_t_complex_gridfun => rset_t_complex_gridfun_tensorial_1D
        procedure :: norm_real_gridfun => norm_real_gridfun_tensorial_1D
        procedure :: norm_complex_gridfun => norm_complex_gridfun_tensorial_1D
!        procedure :: inner_product_real_gridfun => inner_product_real_gridfun_tensorial_1D
!        procedure :: inner_product_complex_gridfun => inner_product_complex_gridfun_tensorial_1D
    end type grid_tensorial_1D
    
!    interface grid_tensorial_1D
!        module procedure new_grid_tensorial_1D 
!    end interface grid_tensorial_1D


    type, extends(grid_2D) :: grid_tensorial_2D
        integer :: nx
        real(kind=prec), allocatable :: nodes_x(:)
        real(kind=prec), allocatable :: weights_x(:)
        integer :: ny
        real(kind=prec), allocatable :: nodes_y(:)
        !real(kind=prec), allocatable :: weights_y(:)
        real(kind=prec), pointer :: weights_y(:)
    contains
        procedure :: set_real_gridfun => set_real_gridfun_tensorial_2D
        procedure :: set_complex_gridfun => set_complex_gridfun_tensorial_2D
        procedure :: rset_complex_gridfun => rset_complex_gridfun_tensorial_2D
        procedure :: set_t_real_gridfun => set_t_real_gridfun_tensorial_2D
        procedure :: set_t_complex_gridfun => set_t_complex_gridfun_tensorial_2D
        procedure :: rset_t_complex_gridfun => rset_t_complex_gridfun_tensorial_2D
        procedure :: norm_real_gridfun => norm_real_gridfun_tensorial_2D
        procedure :: norm_complex_gridfun => norm_complex_gridfun_tensorial_2D
!        procedure :: inner_product_real_gridfun => inner_product_real_gridfun_tensorial_2D
!        procedure :: inner_product_complex_gridfun => inner_product_complex_gridfun_tensorial_2D
    end type grid_tensorial_2D
    
!    interface grid_tensorial_2D
!        module procedure new_grid_tensorial_2D 
!    end interface grid_tensorial_2D


    type, extends(grid_3D) :: grid_tensorial_3D
        integer :: nx
        real(kind=prec), allocatable :: nodes_x(:)
        real(kind=prec), allocatable :: weights_x(:)
        integer :: ny
        real(kind=prec), allocatable :: nodes_y(:)
        real(kind=prec), allocatable :: weights_y(:)
        integer :: nz
        real(kind=prec), allocatable :: nodes_z(:)
        !real(kind=prec), allocatable :: weights_z(:)
        real(kind=prec), pointer :: weights_z(:)
    contains
        procedure :: set_real_gridfun => set_real_gridfun_tensorial_3D
        procedure :: set_complex_gridfun => set_complex_gridfun_tensorial_3D
        procedure :: rset_complex_gridfun => rset_complex_gridfun_tensorial_3D
        procedure :: set_t_real_gridfun => set_t_real_gridfun_tensorial_3D
        procedure :: set_t_complex_gridfun => set_t_complex_gridfun_tensorial_3D
        procedure :: rset_t_complex_gridfun => rset_t_complex_gridfun_tensorial_3D
        procedure :: norm_real_gridfun => norm_real_gridfun_tensorial_3D
        procedure :: norm_complex_gridfun => norm_complex_gridfun_tensorial_3D
!        procedure :: inner_product_real_gridfun => inner_product_real_gridfun_tensorial_3D
!        procedure :: inner_product_complex_gridfun => inner_product_complex_gridfun_tensorial_3D
    end type grid_tensorial_3D
    
!    interface grid_tensorial_3D
!        module procedure new_grid_tensorial_3D 
!    end interface grid_tensorial_3D


    type, extends(grid_2D) :: grid_polar_2D
        integer :: ntheta
        integer :: nr
        integer :: nthetamin, nthetamax
        integer :: nrmin, nrmax
        real(kind=prec) :: dtheta
        real(kind=prec), allocatable :: nodes_r(:)
        !real(kind=prec), allocatable :: weights_r(:)
        real(kind=prec), pointer :: weights_r(:)
        real(kind=prec), allocatable :: nodes_theta(:)
        real(kind=prec), allocatable :: nodes_x(:,:)
        real(kind=prec), allocatable :: nodes_y(:,:)
    contains
        procedure :: compute_nodes_xy => compute_nodes_xy_polar_2D
        procedure :: set_real_gridfun => set_real_gridfun_polar_2D
        procedure :: set_complex_gridfun => set_complex_gridfun_polar_2D
        procedure :: rset_complex_gridfun => rset_complex_gridfun_polar_2D
        procedure :: set_t_real_gridfun => set_t_real_gridfun_polar_2D
        procedure :: set_t_complex_gridfun => set_t_complex_gridfun_polar_2D
        procedure :: rset_t_complex_gridfun => rset_t_complex_gridfun_polar_2D
        procedure :: norm_real_gridfun => norm_real_gridfun_polar_2D
        procedure :: norm_complex_gridfun => norm_complex_gridfun_polar_2D
        procedure :: inner_product_real_gridfun => inner_product_real_gridfun_polar_2D
        procedure :: inner_product_complex_gridfun => inner_product_complex_gridfun_polar_2D
    end type grid_polar_2D
    
!    interface grid_polar_2D
!        module procedure new_grid_polar_2D
!    end interface grid_polar_2D

    type, extends(grid_3D) :: grid_cylindrical_3D
        integer :: ntheta
        integer :: nr
        integer :: nthetamin, nthetamax
        integer :: nrmin, nrmax
        integer :: nzmin, nzmax
        real(kind=prec) :: dtheta
        real(kind=prec), allocatable :: nodes_r(:)
        !real(kind=prec), allocatable :: weights_r(:)
        real(kind=prec), allocatable :: weights_r(:)
        integer :: nz
        real(kind=prec), allocatable :: nodes_z(:)
        !real(kind=prec), allocatable :: weights_z(:)
        real(kind=prec), pointer :: weights_z(:)
        real(kind=prec), allocatable :: nodes_theta(:)
        real(kind=prec), allocatable :: nodes_x(:,:)
        real(kind=prec), allocatable :: nodes_y(:,:)
    contains
        procedure :: compute_nodes_xy => compute_nodes_xy_cylindrical_3D 
        procedure :: set_real_gridfun => set_real_gridfun_cylindrical_3D
        procedure :: set_complex_gridfun => set_complex_gridfun_cylindrical_3D
        procedure :: rset_complex_gridfun => rset_complex_gridfun_cylindrical_3D
        procedure :: set_t_real_gridfun => set_t_real_gridfun_cylindrical_3D
        procedure :: set_t_complex_gridfun => set_t_complex_gridfun_cylindrical_3D
        procedure :: rset_t_complex_gridfun => rset_t_complex_gridfun_cylindrical_3D
        procedure :: norm_real_gridfun => norm_real_gridfun_cylindrical_3D
        procedure :: norm_complex_gridfun => norm_complex_gridfun_cylindrical_3D
        procedure :: inner_product_real_gridfun => inner_product_real_gridfun_cylindrical_3D
        procedure :: inner_product_complex_gridfun => inner_product_complex_gridfun_cylindrical_3D
    end type grid_cylindrical_3D


    
!    interface grid_cylindrical_3D
!        module procedure new_grid_cylindrical_3D
!    end interface grid_cylindrical_3D


contains     


    subroutine allocate_real_gridfun_1D(this, u)
        class(grid_1D) :: this
        real(kind=prec), intent(inout), pointer :: u(:)
        allocate ( u(this%n1min:this%n1max) )
    end subroutine allocate_real_gridfun_1D

    subroutine allocate_complex_gridfun_1D(this, u)
        class(grid_1D) :: this
        complex(kind=prec), intent(inout), pointer :: u(:)
        allocate ( u(this%n1min:this%n1max) )
    end subroutine allocate_complex_gridfun_1D

    subroutine allocate_real_gridfun_2D(this, u)
        class(grid_2D) :: this
        real(kind=prec), intent(inout), pointer :: u(:,:)
        allocate ( u(this%n1min:this%n1max, this%n2min:this%n2max) )
    end subroutine allocate_real_gridfun_2D

    subroutine allocate_complex_gridfun_2D(this, u)
        class(grid_2D) :: this
        complex(kind=prec), intent(inout), pointer :: u(:,:)
        allocate ( u(this%n1min:this%n1max, this%n2min:this%n2max) )
    end subroutine allocate_complex_gridfun_2D

    subroutine allocate_real_gridfun_3D(this, u)
        class(grid_3D) :: this
        real(kind=prec), intent(inout), pointer :: u(:,:,:)
        allocate ( u(this%n1min:this%n1max, this%n2min:this%n2max, this%n3min:this%n3max) )
    end subroutine allocate_real_gridfun_3D

    subroutine allocate_complex_gridfun_3D(this, u)
        class(grid_3D) :: this
        complex(kind=prec), intent(inout), pointer :: u(:,:,:)
        allocate ( u(this%n1min:this%n1max, this%n2min:this%n2max, this%n3min:this%n3max) )
    end subroutine allocate_complex_gridfun_3D


!!!
   subroutine scale_real_gridfun_1D(g, uu, factor)
        class(grid_1D) :: g 
        real(kind=prec), intent(inout), target :: uu(:)
        real(kind=prec), intent(in) :: factor
#ifdef _OPENMP
        real(kind=prec), pointer :: u(:)
        integer :: j
#endif      
#ifndef _OPENMP
        uu = factor*uu
#else
!$OMP PARALLEL DO PRIVATE(j, u) 
         do j=1,n_threads
              u => uu(lbound(uu,1)+g%jj(j-1):lbound(uu,1)+g%jj(j)-1)
              u = factor*u
         end do
!$OMP END PARALLEL DO 
#endif
   end subroutine scale_real_gridfun_1D


   subroutine scale_complex_gridfun_1D(g, uu, factor)
        class(grid_1D) :: g 
        complex(kind=prec), intent(inout), target :: uu(:)
        complex(kind=prec), intent(in) :: factor
#ifdef _OPENMP
        complex(kind=prec), pointer :: u(:)
        integer :: j
#endif      
#ifndef _OPENMP
        if(aimag(factor)==0.0_prec) then
            uu = real(factor, kind=prec)*uu
        else
            uu = factor*uu
        end if
#else
!$OMP PARALLEL DO PRIVATE(j, u) 
         do j=1,n_threads
              u => uu(lbound(uu,1)+g%jj(j-1):lbound(uu,1)+g%jj(j)-1)
              if(aimag(factor)==0.0_prec) then
                  u = real(factor,kind=prec)*u
              else
                  u = factor*u
              end if
         end do
!$OMP END PARALLEL DO 
#endif
   end subroutine scale_complex_gridfun_1D


   subroutine scale_real_gridfun_2D(g, uu, factor)
        class(grid_2D) :: g 
        real(kind=prec), intent(inout), target :: uu(:,:)
        real(kind=prec), intent(in) :: factor
#ifdef _OPENMP
        real(kind=prec), pointer :: u(:,:)
        integer :: j
#endif      
#ifndef _OPENMP
        uu = factor*uu
#else
!$OMP PARALLEL DO PRIVATE(j, u) 
         do j=1,n_threads
              u => uu(:,lbound(uu,2)+g%jj(j-1):lbound(uu,2)+g%jj(j)-1)
              u = factor*u
         end do
!$OMP END PARALLEL DO 
#endif
   end subroutine scale_real_gridfun_2D


   subroutine scale_complex_gridfun_2D(g, uu, factor)
        class(grid_2D) :: g 
        complex(kind=prec), intent(inout), target :: uu(:,:)
        complex(kind=prec), intent(in) :: factor
#ifdef _OPENMP
        complex(kind=prec), pointer :: u(:,:)
        integer :: j
#endif      
#ifndef _OPENMP
        if(aimag(factor)==0.0_prec) then
            uu = real(factor, kind=prec)*uu
        else
            uu = factor*uu
        end if
#else
!$OMP PARALLEL DO PRIVATE(j, u) 
         do j=1,n_threads
              u => uu(:,lbound(uu,2)+g%jj(j-1):lbound(uu,2)+g%jj(j)-1)
              if(aimag(factor)==0.0_prec) then
                  u = real(factor,kind=prec)*u
              else
                  u = factor*u
              end if
         end do
!$OMP END PARALLEL DO 
#endif
   end subroutine scale_complex_gridfun_2D


   subroutine scale_real_gridfun_3D(g, uu, factor)
        class(grid_3D) :: g 
        real(kind=prec), intent(inout), target :: uu(:,:,:)
        real(kind=prec), intent(in) :: factor
#ifdef _OPENMP
        real(kind=prec), pointer :: u(:,:,:)
        integer :: j
#endif      
#ifndef _OPENMP
        uu = factor*uu
#else
!$OMP PARALLEL DO PRIVATE(j, u) 
         do j=1,n_threads
              u => uu(:,:,lbound(uu,3)+g%jj(j-1):lbound(uu,3)+g%jj(j)-1)
              u = factor*u
         end do
!$OMP END PARALLEL DO 
#endif
   end subroutine scale_real_gridfun_3D


   subroutine scale_complex_gridfun_3D(g, uu, factor)
        class(grid_3D) :: g 
        complex(kind=prec), intent(inout), target :: uu(:,:,:)
        complex(kind=prec), intent(in) :: factor
#ifdef _OPENMP
        complex(kind=prec), pointer :: u(:,:,:)
        integer :: j
#endif      
#ifndef _OPENMP
        if(aimag(factor)==0.0_prec) then
            uu = real(factor, kind=prec)*uu
        else
            uu = factor*uu
        end if
#else
!$OMP PARALLEL DO PRIVATE(j, u) 
         do j=1,n_threads
              u => uu(:,:,lbound(uu,3)+g%jj(j-1):lbound(uu,3)+g%jj(j)-1)
              if(aimag(factor)==0.0_prec) then
                  u = real(factor,kind=prec)*u
              else
                  u = factor*u
              end if
         end do
!$OMP END PARALLEL DO 
#endif
   end subroutine scale_complex_gridfun_3D


   subroutine axpy_real_gridfun_1D(g, u_this, u_other, factor)
        class(grid_1D) :: g 
        real(kind=prec), intent(in), target :: u_other(:)
        real(kind=prec), intent(inout), target :: u_this(:)
        real(kind=prec), intent(in) :: factor
#ifdef _OPENMP
        real(kind=prec), pointer :: u1(:)
        real(kind=prec), pointer :: u2(:)
        integer :: j
#endif        
#ifndef _OPENMP
        if(factor==1.0_prec) then
            u_this = u_this + u_other
        elseif(factor==-1.0_prec) then
             u_this = u_this - u_other
        else
            u_this = u_this + factor*u_other
        end if   
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
         do j=1,n_threads
              u1 => u_this(lbound(u_this,1)+g%jj(j-1):lbound(u_this,1)+g%jj(j)-1)
              u2 => u_other(lbound(u_other,1)+g%jj(j-1):lbound(u_other,1)+g%jj(j)-1)
              if(factor==1.0_prec) then
                 u1 = u1 + u2
              elseif(factor==-1.0_prec) then
                 u1 = u1 - u2
              else
                 u1 = u1 + factor*u2
              end if   
         end do
!$OMP END PARALLEL DO 
#endif
   end subroutine axpy_real_gridfun_1D


   subroutine axpy_complex_gridfun_1D(g, u_this, u_other, factor)
        class(grid_1D) :: g 
        complex(kind=prec), intent(in), target :: u_other(:)
        complex(kind=prec), intent(inout), target :: u_this(:)
        complex(kind=prec), intent(in) :: factor
#ifdef _OPENMP
        complex(kind=prec), pointer :: u1(:)
        complex(kind=prec), pointer :: u2(:)
        integer :: j
#endif        
#ifndef _OPENMP
        if(aimag(factor)==0.0_prec) then
            if(real(factor,prec)==1.0_prec) then
               u_this = u_this + u_other
            elseif(real(factor,prec)==-1.0_prec) then
               u_this = u_this - u_other
            else
               u_this = u_this + real(factor, kind=prec)*u_other
            end if   
        else
            u_this = u_this + factor*u_other
        end if
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
         do j=1,n_threads
              u1 => u_this(lbound(u_this,1)+g%jj(j-1):lbound(u_this,1)+g%jj(j)-1)
              u2 => u_other(lbound(u_other,1)+g%jj(j-1):lbound(u_other,1)+g%jj(j)-1)
              if(aimag(factor)==0.0_prec) then
                  if(real(factor,prec)==1.0_prec) then
                     u1 = u1 + u2
                  elseif(real(factor,prec)==-1.0_prec) then
                     u1 = u1 - u2
                  else
                     u1 = u1 + real(factor, kind=prec)*u2
                  end if   
              else
                  u1 = u1 + factor*u2
              end if
         end do
!$OMP END PARALLEL DO 
#endif
   end subroutine axpy_complex_gridfun_1D
   

   subroutine axpy_real_gridfun_2D(g, u_this, u_other, factor)
        class(grid_2D) :: g 
        real(kind=prec), intent(in), target :: u_other(:,:)
        real(kind=prec), intent(inout), target :: u_this(:,:)
        real(kind=prec), intent(in) :: factor
#ifdef _OPENMP
        real(kind=prec), pointer :: u1(:,:)
        real(kind=prec), pointer :: u2(:,:)
        integer :: j
#endif        
#ifndef _OPENMP
        if(factor==1.0_prec) then
            u_this = u_this + u_other
        elseif(factor==-1.0_prec) then
             u_this = u_this - u_other
        else
            u_this = u_this + factor*u_other
        end if   
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
         do j=1,n_threads
              u1 => u_this(:,lbound(u_this,2)+g%jj(j-1):lbound(u_this,2)+g%jj(j)-1)
              u2 => u_other(:,lbound(u_other,2)+g%jj(j-1):lbound(u_other,2)+g%jj(j)-1)
              if(factor==1.0_prec) then
                 u1 = u1 + u2
              elseif(factor==-1.0_prec) then
                 u1 = u1 - u2
              else
                 u1 = u1 + factor*u2
              end if   
         end do
!$OMP END PARALLEL DO 
#endif
   end subroutine axpy_real_gridfun_2D


   subroutine axpy_complex_gridfun_2D(g, u_this, u_other, factor)
        class(grid_2D) :: g 
        complex(kind=prec), intent(in), target :: u_other(:,:)
        complex(kind=prec), intent(inout), target :: u_this(:,:)
        complex(kind=prec), intent(in) :: factor
#ifdef _OPENMP
        complex(kind=prec), pointer :: u1(:,:)
        complex(kind=prec), pointer :: u2(:,:)
        integer :: j
#endif        
#ifndef _OPENMP
        if(aimag(factor)==0.0_prec) then
            if(real(factor,prec)==1.0_prec) then
               u_this = u_this + u_other
            elseif(real(factor,prec)==-1.0_prec) then
               u_this = u_this - u_other
            else
               u_this = u_this + real(factor, kind=prec)*u_other
            end if   
        else
            u_this = u_this + factor*u_other
        end if
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
         do j=1,n_threads
              u1 => u_this(:,lbound(u_this,2)+g%jj(j-1):lbound(u_this,2)+g%jj(j)-1)
              u2 => u_other(:,lbound(u_other,2)+g%jj(j-1):lbound(u_other,2)+g%jj(j)-1)
              if(aimag(factor)==0.0_prec) then
                  if(real(factor,prec)==1.0_prec) then
                     u1 = u1 + u2
                  elseif(real(factor,prec)==-1.0_prec) then
                     u1 = u1 - u2
                  else
                     u1 = u1 + real(factor, kind=prec)*u2
                  end if   
              else
                  u1 = u1 + factor*u2
              end if
         end do
!$OMP END PARALLEL DO 
#endif
   end subroutine axpy_complex_gridfun_2D


   subroutine axpy_real_gridfun_3D(g, u_this, u_other, factor)
        class(grid_3D) :: g 
        real(kind=prec), intent(in), target :: u_other(:,:,:)
        real(kind=prec), intent(inout), target :: u_this(:,:,:)
        real(kind=prec), intent(in) :: factor
#ifdef _OPENMP
        real(kind=prec), pointer :: u1(:,:,:)
        real(kind=prec), pointer :: u2(:,:,:)
        integer :: j
#endif        
#ifndef _OPENMP
        if(factor==1.0_prec) then
            u_this = u_this + u_other
        elseif(factor==-1.0_prec) then
             u_this = u_this - u_other
        else
            u_this = u_this + factor*u_other
        end if   
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
         do j=1,n_threads
              u1 => u_this(:,:,lbound(u_this,3)+g%jj(j-1):lbound(u_this,3)+g%jj(j)-1)
              u2 => u_other(:,:,lbound(u_other,3)+g%jj(j-1):lbound(u_other,3)+g%jj(j)-1)
              if(factor==1.0_prec) then
                 u1 = u1 + u2
              elseif(factor==-1.0_prec) then
                 u1 = u1 - u2
              else
                 u1 = u1 + factor*u2
              end if   
         end do
!$OMP END PARALLEL DO 
#endif
   end subroutine axpy_real_gridfun_3D


   subroutine axpy_complex_gridfun_3D(g, u_this, u_other, factor)
        class(grid_3D) :: g 
        complex(kind=prec), intent(in), target :: u_other(:,:,:)
        complex(kind=prec), intent(inout), target :: u_this(:,:,:)
        complex(kind=prec), intent(in) :: factor
#ifdef _OPENMP
        complex(kind=prec), pointer :: u1(:,:,:)
        complex(kind=prec), pointer :: u2(:,:,:)
        integer :: j
#endif        
#ifndef _OPENMP
        if(aimag(factor)==0.0_prec) then
            if(real(factor,prec)==1.0_prec) then
               u_this = u_this + u_other
            elseif(real(factor,prec)==-1.0_prec) then
               u_this = u_this - u_other
            else
               u_this = u_this + real(factor, kind=prec)*u_other
            end if   
        else
            u_this = u_this + factor*u_other
        end if
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
         do j=1,n_threads
              u1 => u_this(:,:,lbound(u_this,3)+g%jj(j-1):lbound(u_this,3)+g%jj(j)-1)
              u2 => u_other(:,:,lbound(u_other,3)+g%jj(j-1):lbound(u_other,3)+g%jj(j)-1)
              if(aimag(factor)==0.0_prec) then
                  if(real(factor,prec)==1.0_prec) then
                     u1 = u1 + u2
                  elseif(real(factor,prec)==-1.0_prec) then
                     u1 = u1 - u2
                  else
                     u1 = u1 + real(factor, kind=prec)*u2
                  end if   
              else
                  u1 = u1 + factor*u2
              end if
         end do
!$OMP END PARALLEL DO 
#endif
   end subroutine axpy_complex_gridfun_3D

!!!

   subroutine copy_real_gridfun_1D(g, u_target, u_source)
        class(grid_1D) :: g 
        real(kind=prec), intent(in), target :: u_source(:)
        real(kind=prec), intent(out), target :: u_target(:)
#ifdef _OPENMP
        real(kind=prec), pointer :: u1(:)
        real(kind=prec), pointer :: u2(:)
        integer :: j
#endif        
#ifndef _OPENMP
         u_target = u_source
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
         do j=1,n_threads
              u1 => u_target(lbound(u_target,1)+g%jj(j-1):lbound(u_target,1)+g%jj(j)-1)
              u2 => u_source(lbound(u_source,1)+g%jj(j-1):lbound(u_source,1)+g%jj(j)-1)
              u1 = u2
         end do
!$OMP END PARALLEL DO 
#endif
   end subroutine copy_real_gridfun_1D


   subroutine copy_complex_gridfun_1D(g, u_target, u_source)
        class(grid_1D) :: g 
        complex(kind=prec), intent(in), target :: u_source(:)
        complex(kind=prec), intent(out), target :: u_target(:)
#ifdef _OPENMP
        complex(kind=prec), pointer :: u1(:)
        complex(kind=prec), pointer :: u2(:)
        integer :: j
#endif        
#ifndef _OPENMP
         u_target = u_source
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
         do j=1,n_threads
              u1 => u_target(lbound(u_target,1)+g%jj(j-1):lbound(u_target,1)+g%jj(j)-1)
              u2 => u_source(lbound(u_source,1)+g%jj(j-1):lbound(u_source,1)+g%jj(j)-1)
              u1 = u2
         end do
!$OMP END PARALLEL DO 
#endif
   end subroutine copy_complex_gridfun_1D


   subroutine copy_real_gridfun_2D(g, u_target, u_source)
        class(grid_2D) :: g 
        real(kind=prec), intent(in), target :: u_source(:,:)
        real(kind=prec), intent(out), target :: u_target(:,:)
#ifdef _OPENMP
        real(kind=prec), pointer :: u1(:,:)
        real(kind=prec), pointer :: u2(:,:)
        integer :: j
#endif        
#ifndef _OPENMP
         u_target = u_source
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
         do j=1,n_threads
              u1 => u_target(:,lbound(u_target,2)+g%jj(j-1):lbound(u_target,2)+g%jj(j)-1)
              u2 => u_source(:,lbound(u_source,2)+g%jj(j-1):lbound(u_source,2)+g%jj(j)-1)
              u1 = u2
         end do
!$OMP END PARALLEL DO 
#endif
   end subroutine copy_real_gridfun_2D


   subroutine copy_complex_gridfun_2D(g, u_target, u_source)
        class(grid_2D) :: g 
        complex(kind=prec), intent(in), target :: u_source(:,:)
        complex(kind=prec), intent(out), target :: u_target(:,:)
#ifdef _OPENMP
        complex(kind=prec), pointer :: u1(:,:)
        complex(kind=prec), pointer :: u2(:,:)
        integer :: j
#endif        
#ifndef _OPENMP
         u_target = u_source
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
         do j=1,n_threads
              u1 => u_target(:,lbound(u_target,2)+g%jj(j-1):lbound(u_target,2)+g%jj(j)-1)
              u2 => u_source(:,lbound(u_source,2)+g%jj(j-1):lbound(u_source,2)+g%jj(j)-1)
              u1 = u2
         end do
!$OMP END PARALLEL DO 
#endif
   end subroutine copy_complex_gridfun_2D


   subroutine copy_real_gridfun_3D(g, u_target, u_source)
        class(grid_3D) :: g 
        real(kind=prec), intent(in), target :: u_source(:,:,:)
        real(kind=prec), intent(out), target :: u_target(:,:,:)
#ifdef _OPENMP
        real(kind=prec), pointer :: u1(:,:,:)
        real(kind=prec), pointer :: u2(:,:,:)
        integer :: j
#endif        
#ifndef _OPENMP
         u_target = u_source
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
         do j=1,n_threads
              u1 => u_target(:,:,lbound(u_target,3)+g%jj(j-1):lbound(u_target,3)+g%jj(j)-1)
              u2 => u_source(:,:,lbound(u_source,3)+g%jj(j-1):lbound(u_source,3)+g%jj(j)-1)
              u1 = u2
         end do
!$OMP END PARALLEL DO 
#endif
   end subroutine copy_real_gridfun_3D


   subroutine copy_complex_gridfun_3D(g, u_target, u_source)
        class(grid_3D) :: g 
        complex(kind=prec), intent(in), target :: u_source(:,:,:)
        complex(kind=prec), intent(out), target :: u_target(:,:,:)
#ifdef _OPENMP
        complex(kind=prec), pointer :: u1(:,:,:)
        complex(kind=prec), pointer :: u2(:,:,:)
        integer :: j
#endif        
#ifndef _OPENMP
         u_target = u_source
#else
!$OMP PARALLEL DO PRIVATE(j, u1, u2) 
         do j=1,n_threads
              u1 => u_target(:,:,lbound(u_target,3)+g%jj(j-1):lbound(u_target,3)+g%jj(j)-1)
              u2 => u_source(:,:,lbound(u_source,3)+g%jj(j-1):lbound(u_source,3)+g%jj(j)-1)
              u1 = u2
         end do
!$OMP END PARALLEL DO 
#endif
   end subroutine copy_complex_gridfun_3D




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   
   subroutine set_real_gridfun_equidistant_1D(this, u, f)
        class(grid_equidistant_1D) :: this
        real(kind=prec), intent(inout) :: u(this%n1min:this%n1max)
        real(kind=prec), external :: f

        integer :: ix
        real(kind=prec) :: x

!$OMP PARALLEL DO PRIVATE(ix, x)
        do ix = this%n1min,this%n1max
            x = this%xmin + this%dx*ix  
            u(ix) = f(x)
        end do
!$OMP END PARALLEL DO
    end subroutine set_real_gridfun_equidistant_1D

    subroutine set_complex_gridfun_equidistant_1D(this, u, f)
        class(grid_equidistant_1D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max)
        complex(kind=prec), external :: f

        integer :: ix
        real(kind=prec) :: x
!TODO: use nodes_x instead of computing x ...
!$OMP PARALLEL DO PRIVATE(ix, x)
        do ix = this%n1min,this%n1max
            x = this%xmin + this%dx*ix    
            u(ix) = f(x)
        end do
!$OMP END PARALLEL DO
    end subroutine set_complex_gridfun_equidistant_1D

    subroutine rset_complex_gridfun_equidistant_1D(this, u, f)
        class(grid_equidistant_1D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max)
        real(kind=prec), external :: f

        integer :: ix
        real(kind=prec) :: x

!$OMP PARALLEL DO PRIVATE(ix, x)
        do ix = this%n1min,this%n1max
            x = this%xmin + this%dx*ix    
            u(ix) = cmplx(f(x), kind=prec)
        end do
!$OMP END PARALLEL DO
    end subroutine rset_complex_gridfun_equidistant_1D

   subroutine set_t_real_gridfun_equidistant_1D(this, u, f, t)
        class(grid_equidistant_1D) :: this
        real(kind=prec), intent(inout) :: u(this%n1min:this%n1max)
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t

        integer :: ix
        real(kind=prec) :: x

!$OMP PARALLEL DO PRIVATE(ix, x)
        do ix = this%n1min,this%n1max
            x = this%xmin + this%dx*ix  
            u(ix) = f(x, t)
        end do
!$OMP END PARALLEL DO
    end subroutine set_t_real_gridfun_equidistant_1D

    subroutine set_t_complex_gridfun_equidistant_1D(this, u, f, t)
        class(grid_equidistant_1D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max)
        complex(kind=prec), external :: f
        real(kind=prec), intent(in) :: t

        integer :: ix
        real(kind=prec) :: x
!TODO: use nodes_x instead of computing x ...
!$OMP PARALLEL DO PRIVATE(ix, x)
        do ix = this%n1min,this%n1max
            x = this%xmin + this%dx*ix    
            u(ix) = f(x, t)
        end do
!$OMP END PARALLEL DO
    end subroutine set_t_complex_gridfun_equidistant_1D

    subroutine rset_t_complex_gridfun_equidistant_1D(this, u, f, t)
        class(grid_equidistant_1D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max)
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t

        integer :: ix
        real(kind=prec) :: x

!$OMP PARALLEL DO PRIVATE(ix, x)
        do ix = this%n1min,this%n1max
            x = this%xmin + this%dx*ix    
            u(ix) = cmplx(f(x, t), kind=prec)
        end do
!$OMP END PARALLEL DO
    end subroutine rset_t_complex_gridfun_equidistant_1D


    function norm_real_gridfun_equidistant_1D(this, u) result(n)
        class(grid_equidistant_1D) :: this
        !real(kind=prec), intent(in) :: u(this%n1min:this%n1max)
        real(kind=prec), intent(in), target :: u(:)
        real(kind=prec) :: n
        real(kind=prec) :: s
#ifdef _OPENMP
        real(kind=prec), pointer :: v(:)
        integer :: j
#endif         
#ifdef _MPI_
        real(kind=prec) :: s1
        integer :: ierr
#endif         
#ifndef _OPENMP
        s = sum(u**2)
#else
        s = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v) REDUCTION(+:s)
        do j=1,n_threads
            v => u(lbound(u,1)+this%jj(j-1):lbound(u,1)+this%jj(j)-1)        
            s = s + sum(v**2)
        end do    
!$OMP END PARALLEL DO 
#endif
#ifdef _MPI_
        call MPI_Reduce(s, s1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        s = s1
        call MPI_Bcast(s, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#endif
        n = sqrt(s*this%dx)
    end function norm_real_gridfun_equidistant_1D


    function inner_product_real_gridfun_equidistant_1D(this, u1, u2) result(n)
        class(grid_equidistant_1D) :: this
        !real(kind=prec), intent(in) :: u(this%n1min:this%n1max)
        real(kind=prec), intent(in), target :: u1(:)
        real(kind=prec), intent(in), target :: u2(:)
        real(kind=prec) :: n
        real(kind=prec) :: s
#ifdef _OPENMP
        real(kind=prec), pointer :: v1(:)
        real(kind=prec), pointer :: v2(:)
        integer :: j
#endif         
#ifdef _MPI_
        real(kind=prec) :: s1
        integer :: ierr
#endif         
#ifndef _OPENMP
        s = sum(u1*u2)
#else
        s = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v1, v2) REDUCTION(+:s)
        do j=1,n_threads
            v1 => u1(lbound(u1,1)+this%jj(j-1):lbound(u1,1)+this%jj(j)-1)        
            v2 => u2(lbound(u2,1)+this%jj(j-1):lbound(u2,1)+this%jj(j)-1)        
            s = s + sum(v1*v2)
        end do    
!$OMP END PARALLEL DO 
#endif
#ifdef _MPI_
        call MPI_Reduce(s, s1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        s = s1
        call MPI_Bcast(s, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#endif
        n = s*this%dx
    end function inner_product_real_gridfun_equidistant_1D



    function norm_complex_gridfun_equidistant_1D(this, u) result(n)
        class(grid_equidistant_1D) :: this
        !complex(kind=prec), intent(in) :: u(this%n1min:this%n1max)
        complex(kind=prec), intent(in), target :: u(:)
        real(kind=prec) :: n
        real(kind=prec) :: s
#ifdef _OPENMP
        complex(kind=prec), pointer :: v(:)
        integer :: j
#endif              
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) ::  s1
#endif        
#ifndef _OPENMP
        s = sum(real(u,kind=prec)**2 +aimag(u)**2)
#else
        s = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v) REDUCTION(+:s)
        do j=1,n_threads
            v => u(lbound(u,1)+this%jj(j-1):lbound(u,1)+this%jj(j)-1)        
            s = s + sum(real(v,kind=prec)**2 +aimag(v)**2)
        end do    
!$OMP END PARALLEL DO 
#endif
#ifdef _MPI_
        call MPI_Reduce(s, s1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        s = s1
        call MPI_Bcast(s, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#endif        
        n = sqrt(s*this%dx)
    end function norm_complex_gridfun_equidistant_1D


    function inner_product_complex_gridfun_equidistant_1D(this, u1, u2) result(n)
        class(grid_equidistant_1D) :: this
        !complex(kind=prec), intent(in) :: u(this%n1min:this%n1max)
        complex(kind=prec), intent(in), target :: u1(:)
        complex(kind=prec), intent(in), target :: u2(:)
        complex(kind=prec) :: n
        complex(kind=prec) :: s
#ifdef _OPENMP
        complex(kind=prec), pointer :: v1(:)
        complex(kind=prec), pointer :: v2(:)
        integer :: j
#endif              
#ifdef _MPI_
        integer :: ierr
        complex(kind=prec) ::  s1
#endif        
#ifndef _OPENMP
        s = sum(conjg(u1)*u2)
#else
        s = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v1, v2) REDUCTION(+:s)
        do j=1,n_threads
            v1 => u1(lbound(u1,1)+this%jj(j-1):lbound(u1,1)+this%jj(j)-1)        
            v2 => u2(lbound(u2,1)+this%jj(j-1):lbound(u2,1)+this%jj(j)-1)        
            s = s + sum(conjg(v1)*v2)
        end do    
!$OMP END PARALLEL DO 
#endif
#ifdef _MPI_
        call MPI_Reduce(s, s1, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        s = s1
        call MPI_Bcast(s, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr) 
#endif        
        n = s*this%dx
    end function inner_product_complex_gridfun_equidistant_1D


    subroutine set_real_gridfun_equidistant_2D(this, u, f)
        class(grid_equidistant_2D) :: this
        real(kind=prec), intent(out) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
        real(kind=prec), external :: f

        integer :: ix, iy
        real(kind=prec) :: x, y

!$OMP PARALLEL DO PRIVATE(ix, x, iy, y)
        do iy = this%n2min, this%n2max
            y = this%ymin + this%dy*iy    
            do ix = this%n1min,this%n1max
                x = this%xmin + this%dx*ix   
                u(ix, iy) = f(x, y)
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine set_real_gridfun_equidistant_2D

    subroutine set_complex_gridfun_equidistant_2D(this, u, f)
        class(grid_equidistant_2D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
        complex(kind=prec), external :: f

        integer :: ix, iy
        real(kind=prec) :: x, y

!$OMP PARALLEL DO PRIVATE(ix, x, iy, y)
        do iy = this%n2min, this%n2max
            y = this%ymin + this%dy*iy    
            do ix = this%n1min,this%n1max
                x = this%xmin + this%dx*ix    
                u(ix, iy) = f(x, y)
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine set_complex_gridfun_equidistant_2D

    subroutine rset_complex_gridfun_equidistant_2D(this, u, f)
        class(grid_equidistant_2D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
        real(kind=prec), external :: f

        integer :: ix, iy
        real(kind=prec) :: x, y
!$OMP PARALLEL DO PRIVATE(ix, x, iy, y)
        do iy = this%n2min, this%n2max
            y = this%ymin + this%dy*iy    
            do ix = this%n1min,this%n1max
                x = this%xmin + this%dx*ix    
                u(ix, iy) = cmplx(f(x, y), kind=prec)
            end do
        end do
!$OMP END PARALLEL DO
end subroutine rset_complex_gridfun_equidistant_2D

    subroutine set_t_real_gridfun_equidistant_2D(this, u, f, t)
        class(grid_equidistant_2D) :: this
        real(kind=prec), intent(out) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t

        integer :: ix, iy
        real(kind=prec) :: x, y

!$OMP PARALLEL DO PRIVATE(ix, x, iy, y)
        do iy = this%n2min, this%n2max
            y = this%ymin + this%dy*iy    
            do ix = this%n1min,this%n1max
                x = this%xmin + this%dx*ix   
                u(ix, iy) = f(x, y, t)
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine set_t_real_gridfun_equidistant_2D

    subroutine set_t_complex_gridfun_equidistant_2D(this, u, f, t)
        class(grid_equidistant_2D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
        complex(kind=prec), external :: f
        real(kind=prec), intent(in) :: t

        integer :: ix, iy
        real(kind=prec) :: x, y

!$OMP PARALLEL DO PRIVATE(ix, x, iy, y)
        do iy = this%n2min, this%n2max
            y = this%ymin + this%dy*iy    
            do ix = this%n1min,this%n1max
                x = this%xmin + this%dx*ix    
                u(ix, iy) = f(x, y, t)
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine set_t_complex_gridfun_equidistant_2D

    subroutine rset_t_complex_gridfun_equidistant_2D(this, u, f, t)
        class(grid_equidistant_2D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t

        integer :: ix, iy
        real(kind=prec) :: x, y
!$OMP PARALLEL DO PRIVATE(ix, x, iy, y)
        do iy = this%n2min, this%n2max
            y = this%ymin + this%dy*iy    
            do ix = this%n1min,this%n1max
                x = this%xmin + this%dx*ix    
                u(ix, iy) = cmplx(f(x, y, t), kind=prec)
            end do
        end do
!$OMP END PARALLEL DO
end subroutine rset_t_complex_gridfun_equidistant_2D


function norm_real_gridfun_equidistant_2D(this, u) result(n)
class(grid_equidistant_2D) :: this
!real(kind=prec), intent(in) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
real(kind=prec), intent(in), target :: u(:,:)
real(kind=prec) :: n
real(kind=prec) :: s
#ifdef _OPENMP
real(kind=prec), pointer :: v(:,:)
integer :: j
#endif         
#ifdef _MPI_
real(kind=prec) :: s1
integer :: ierr
#endif         
#ifndef _OPENMP
s = sum(u**2)
#else
s = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v) REDUCTION(+:s)
do j=1,n_threads
    v => u(:,lbound(u,2)+this%jj(j-1):lbound(u,2)+this%jj(j)-1)        
    s = s + sum(v**2)
end do    
!$OMP END PARALLEL DO 
#endif
#ifdef _MPI_
call MPI_Reduce(s, s1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
s = s1
call MPI_Bcast(s, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#endif
n = sqrt(s*this%dx*this%dy)
end function norm_real_gridfun_equidistant_2D


function norm_complex_gridfun_equidistant_2D(this, u) result(n)
class(grid_equidistant_2D) :: this
!complex(kind=prec), intent(in) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
complex(kind=prec), intent(in), target :: u(:,:)
real(kind=prec) :: n
real(kind=prec) :: s
#ifdef _OPENMP
complex(kind=prec), pointer :: v(:,: )
integer :: j
#endif              
#ifdef _MPI_
integer :: ierr
real(kind=prec) ::  s1
#endif        
#ifndef _OPENMP
s = sum(real(u,kind=prec)**2 +aimag(u)**2)
#else
s = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v) REDUCTION(+:s)
do j=1,n_threads
    v => u(:,lbound(u,2)+this%jj(j-1):lbound(u,2)+this%jj(j)-1)        
    s = s + sum(real(v,kind=prec)**2 +aimag(v)**2)
end do    
!$OMP END PARALLEL DO 
#endif
#ifdef _MPI_
call MPI_Reduce(s, s1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
s = s1
call MPI_Bcast(s, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#endif        
n = sqrt(s*this%dx*this%dy)
end function norm_complex_gridfun_equidistant_2D


subroutine set_real_gridfun_equidistant_3D(this, u, f)
        class(grid_equidistant_3D) :: this
        real(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        real(kind=prec), external :: f

        integer :: ix, iy, iz
        real(kind=prec) :: x, y, z

!$OMP PARALLEL DO PRIVATE(ix, x, iy, y, iz, z) 
        do iz = this%n3min, this%n3max
            z = this%zmin + this%dz*iz    
            do iy = this%n2min, this%n2max
                y = this%ymin + this%dy*iy    
                do ix = this%n1min,this%n1max
                    x = this%xmin + this%dx*ix    
                    u(ix, iy, iz) = f(x, y, z)
                end do
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine set_real_gridfun_equidistant_3D


    subroutine set_complex_gridfun_equidistant_3D(this, u, f)
        class(grid_equidistant_3D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        complex(kind=prec), external :: f

        integer :: ix, iy, iz
        real(kind=prec) :: x, y, z

!$OMP PARALLEL DO PRIVATE(ix, x, iy, y, iz, z)
        do iz = this%n3min, this%n3max
            z = this%zmin + this%dz*iz    
            do iy = this%n2min, this%n2max
                y = this%ymin + this%dy*iy    
                do ix = this%n1min,this%n1max
                    x = this%xmin + this%dx*ix    
                    u(ix, iy, iz) = f(x, y, z)
                end do
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine set_complex_gridfun_equidistant_3D

    subroutine rset_complex_gridfun_equidistant_3D(this, u, f)
        class(grid_equidistant_3D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        real(kind=prec), external :: f

        integer :: ix, iy, iz
        real(kind=prec) :: x, y, z

!$OMP PARALLEL DO PRIVATE(ix, x, iy, y, iz, z)
        do iz = this%n3min, this%n3max
            z = this%zmin + this%dz*iz    
            do iy = this%n2min, this%n2max
                y = this%ymin + this%dy*iy    
                do ix = this%n1min,this%n1max
                    x = this%xmin + this%dx*ix    
                    u(ix, iy, iz) = cmplx(f(x, y, z), kind=prec)
                end do
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine rset_complex_gridfun_equidistant_3D


    subroutine set_t_real_gridfun_equidistant_3D(this, u, f, t)
        class(grid_equidistant_3D) :: this
        real(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t

        integer :: ix, iy, iz
        real(kind=prec) :: x, y, z

!$OMP PARALLEL DO PRIVATE(ix, x, iy, y, iz, z) 
        do iz = this%n3min, this%n3max
            z = this%zmin + this%dz*iz    
            do iy = this%n2min, this%n2max
                y = this%ymin + this%dy*iy    
                do ix = this%n1min,this%n1max
                    x = this%xmin + this%dx*ix    
                    u(ix, iy, iz) = f(x, y, z, t)
                end do
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine set_t_real_gridfun_equidistant_3D


    subroutine set_t_complex_gridfun_equidistant_3D(this, u, f, t)
        class(grid_equidistant_3D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        complex(kind=prec), external :: f
        real(kind=prec), intent(in) :: t

        integer :: ix, iy, iz
        real(kind=prec) :: x, y, z

!$OMP PARALLEL DO PRIVATE(ix, x, iy, y, iz, z)
        do iz = this%n3min, this%n3max
            z = this%zmin + this%dz*iz    
            do iy = this%n2min, this%n2max
                y = this%ymin + this%dy*iy    
                do ix = this%n1min,this%n1max
                    x = this%xmin + this%dx*ix    
                    u(ix, iy, iz) = f(x, y, z, t)
                end do
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine set_t_complex_gridfun_equidistant_3D

    subroutine rset_t_complex_gridfun_equidistant_3D(this, u, f, t)
        class(grid_equidistant_3D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t

        integer :: ix, iy, iz
        real(kind=prec) :: x, y, z

!$OMP PARALLEL DO PRIVATE(ix, x, iy, y, iz, z)
        do iz = this%n3min, this%n3max
            z = this%zmin + this%dz*iz    
            do iy = this%n2min, this%n2max
                y = this%ymin + this%dy*iy    
                do ix = this%n1min,this%n1max
                    x = this%xmin + this%dx*ix    
                    u(ix, iy, iz) = cmplx(f(x, y, z, t), kind=prec)
                end do
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine rset_t_complex_gridfun_equidistant_3D

    
    function norm_real_gridfun_equidistant_3D(this, u) result(n)
        class(grid_equidistant_3D) :: this
        !real(kind=prec), intent(in) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        real(kind=prec), intent(in), target :: u(:,:,:)
        real(kind=prec) :: n
        real(kind=prec) :: s
#ifdef _OPENMP
        real(kind=prec), pointer :: v(:,:,:)
        integer :: j
#endif         
#ifdef _MPI_
        real(kind=prec) :: s1
        integer :: ierr
#endif         
#ifndef _OPENMP
        s = sum(u**2)
#else
        s = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v) REDUCTION(+:s)
        do j=1,n_threads
            v => u(:,:,lbound(u,3)+this%jj(j-1):lbound(u,3)+this%jj(j)-1)        
            s = s + sum(v**2)
        end do    
!$OMP END PARALLEL DO 
#endif
#ifdef _MPI_
        call MPI_Reduce(s, s1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        s = s1
        call MPI_Bcast(s, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#endif
        n = sqrt(s*this%dx*this%dy*this%dz)
    end function norm_real_gridfun_equidistant_3D


    function norm_complex_gridfun_equidistant_3D(this, u) result(n)
        class(grid_equidistant_3D) :: this
        !complex(kind=prec), intent(in) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        complex(kind=prec), intent(in), target :: u(:,:,:)
        real(kind=prec) :: n
        real(kind=prec) :: s
#ifdef _OPENMP
        complex(kind=prec), pointer :: v(:,:,:)
        integer :: j
#endif              
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) ::  s1
#endif        
#ifndef _OPENMP
        s = sum(real(u,kind=prec)**2 +aimag(u)**2)
#else
        s = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v) REDUCTION(+:s)
        do j=1,n_threads
            v => u(:,:,lbound(u,3)+this%jj(j-1):lbound(u,3)+this%jj(j)-1)        
            s = s + sum(real(v,kind=prec)**2 +aimag(v)**2)
        end do    
!$OMP END PARALLEL DO 
#endif
#ifdef _MPI_
        call MPI_Reduce(s, s1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        s = s1
        call MPI_Bcast(s, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#endif        
        n = sqrt(s*this%dx*this%dy*this%dz)
    end function norm_complex_gridfun_equidistant_3D


    function new_grid_equidistant_1D(nx, xmin, xmax) result(this)
        type(grid_equidistant_1D) :: this
        integer, intent(in) :: nx
        real(kind=prec), intent(in) :: xmin, xmax

        this%nx =nx
        this%xmin = xmin
        this%xmax = xmax

        this%nn1min = 1
        this%nn1max = this%nx         
        this%n1min = this%nn1min
        this%n1max = this%nn1max
        this%m1min = this%nn1min
        this%m1max = this%nn1max

        this%dx = (this%xmax-this%xmin)/this%nx
    end function new_grid_equidistant_1D

    function new_grid_equidistant_2D(nx, xmin, xmax, ny, ymin, ymax) result(this)
        type(grid_equidistant_2D):: this
        integer, intent(in) :: nx, ny
        real(kind=prec), intent(in) :: xmin, xmax, ymin, ymax

        this%nx =nx
        this%xmin = xmin
        this%xmax = xmax
        this%ny =ny
        this%ymin = ymin
        this%ymax = ymax

        this%nn1min = 1
        this%nn1max = this%nx         
        this%n1min = this%nn1min
        this%n1max = this%nn1max
        this%m1min = this%nn1min
        this%m1max = this%nn1max
        this%nn2min = 1
        this%nn2max = this%ny         
        this%n2min = this%nn2min
        this%n2max = this%nn2max
        this%m2min = this%nn2min
        this%m2max = this%nn2max

        this%dx = (this%xmax-this%xmin)/this%nx
        this%dy = (this%ymax-this%ymin)/this%ny
    end function new_grid_equidistant_2D

    function new_grid_equidistant_3D(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax) result(this)
        type(grid_equidistant_3D):: this
        integer, intent(in) :: nx, ny, nz
        real(kind=prec), intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax

        this%nx =nx
        this%xmin = xmin
        this%xmax = xmax
        this%ny =ny
        this%ymin = ymin
        this%ymax = ymax
        this%nz =nz
        this%zmin = zmin
        this%zmax = zmax

        this%nn1min = 1
        this%nn1max = this%nx         
        this%n1min = this%nn1min
        this%n1max = this%nn1max
        this%m1min = this%nn1min
        this%m1max = this%nn1max
        this%nn2min = 1
        this%nn2max = this%ny         
        this%n2min = this%nn2min
        this%n2max = this%nn2max
        this%m2min = this%nn2min
        this%m2max = this%nn2max
        this%nn3min = 1
        this%nn3max = this%nz         
        this%n3min = this%nn3min
        this%n3max = this%nn3max
        this%m3min = this%nn3min
        this%m3max = this%nn3max

        this%dx = (this%xmax-this%xmin)/this%nx
        this%dy = (this%ymax-this%ymin)/this%ny
        this%dz = (this%zmax-this%zmin)/this%nz
    end function new_grid_equidistant_3D


    function norm_real_gridfun_tensorial_1D(this, u) result(n)
        class(grid_tensorial_1D) :: this
        !real(kind=prec), intent(in) :: u(this%n1min:this%n1max)
        real(kind=prec), intent(in), target :: u(:)
        real(kind=prec) :: n
#ifdef _OPENMP
        real(kind=prec), pointer :: v(:)
        real(kind=prec), pointer :: w(:)
        integer :: j
#endif            
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: n1
#endif       
#ifndef _OPENMP
        n = sum(this%weights_x*u**2)
#else
        n = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v, w) REDUCTION(+:n)
        do j=1,n_threads
            v => u(lbound(u,1)+this%jj(j-1):lbound(u,1)+this%jj(j)-1)        
            w => this%weights_x(lbound(this%weights_x,1)+this%jj(j-1):lbound(this%weights_x,1)+this%jj(j)-1)
            n = n + sum(w * v**2)
        end do    
!$OMP END PARALLEL DO 
#endif 
#ifdef _MPI_
        call MPI_Reduce(n, n1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            n = sqrt(n1)
        end if    
        call MPI_Bcast(n, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#else
        n = sqrt(n)
#endif        
    end function norm_real_gridfun_tensorial_1D


    function norm_complex_gridfun_tensorial_1D(this, u) result(n)
        class(grid_tensorial_1D) :: this
        !complex(kind=prec), intent(in) :: u(this%n1min:this%n1max)
        complex(kind=prec), intent(in), target :: u(:)
        real(kind=prec) :: n
#ifdef _OPENMP
        complex(kind=prec), pointer :: v(:)
        real(kind=prec), pointer :: w(:)
        integer :: j
#endif            
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: n1
#endif        
#ifndef _OPENMP
        n = sum(this%weights_x*(real(u,kind=prec)**2 +aimag(u)**2))
#else
        n = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v, w) REDUCTION(+:n)
        do j=1,n_threads
            v => u(lbound(u,1)+this%jj(j-1):lbound(u,1)+this%jj(j)-1)        
            w => this%weights_x(lbound(this%weights_x,1)+this%jj(j-1):lbound(this%weights_x,1)+this%jj(j)-1)
            n = n + sum( w *(real(v,kind=prec)**2 +aimag(v)**2))
        end do    
!$OMP END PARALLEL DO 
#endif 
#ifdef _MPI_
        call MPI_Reduce(n, n1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            n = sqrt(n1)
        end if    
        call MPI_Bcast(n, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#else
        n = sqrt(n)
#endif        
    end function norm_complex_gridfun_tensorial_1D


   function norm_real_gridfun_tensorial_2D(this, u) result(n)
        class(grid_tensorial_2D) :: this
        real(kind=prec), intent(in), target :: u(:,:)
        real(kind=prec) :: n
#ifdef _OPENMP
        real(kind=prec), pointer :: v(:,:)
        real(kind=prec), pointer :: w(:)
        integer :: j
#endif                  
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: n1
#endif        
#ifndef _OPENMP
        n = sum(spread(this%weights_x, 2, this%n2max-this%n2min+1) &
               *spread(this%weights_y, 1, this%n1max-this%n1min+1) &
               *u**2)
#else
        n = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v, w) REDUCTION(+:n)
        do j=1,n_threads
            v => u(:,lbound(u,2)+this%jj(j-1):lbound(u,2)+this%jj(j)-1)        
            w => this%weights_y(lbound(this%weights_y,1)+this%jj(j-1):lbound(this%weights_y,1)+this%jj(j)-1)
            n = n + sum(spread(this%weights_x, 2, this%jj(j)-this%jj(j-1) ) &
                       *spread(w, 1, this%n1max-this%n1min+1) &
                       *v**2)
        end do    
!$OMP END PARALLEL DO 
#endif 
#ifdef _MPI_
        call MPI_Reduce(n, n1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            n = sqrt(n1)
        end if    
        call MPI_Bcast(n, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#else
        n = sqrt(n)
#endif        
    end function norm_real_gridfun_tensorial_2D


    function norm_complex_gridfun_tensorial_2D(this, u) result(n)
        class(grid_tensorial_2D) :: this
        complex(kind=prec), intent(in), target :: u(:,:)
        real(kind=prec) :: n
#ifdef _OPENMP
        complex(kind=prec), pointer :: v(:,:)
        real(kind=prec), pointer :: w(:)
        integer :: j
#endif                    
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: n1
#endif  
#ifndef _OPENMP
        n = sum(spread(this%weights_x, 2, this%n2max-this%n2min+1) &
               *spread(this%weights_y, 1, this%n1max-this%n1min+1) &
               *(real(u,kind=prec)**2 +aimag(u)**2))
#else
        n = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v, w) REDUCTION(+:n)
        do j=1,n_threads
            v => u(:,lbound(u,2)+this%jj(j-1):lbound(u,2)+this%jj(j)-1)        
            w => this%weights_y(lbound(this%weights_y,1)+this%jj(j-1):lbound(this%weights_y,1)+this%jj(j)-1)
            n = n + sum(spread(this%weights_x, 2, this%jj(j)-this%jj(j-1) ) &
                       *spread(w, 1, this%n1max-this%n1min+1) &
                       *(real(v,kind=prec)**2 +aimag(v)**2))
        end do    
!$OMP END PARALLEL DO 
#endif
#ifdef _MPI_
        call MPI_Reduce(n, n1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            n = sqrt(n1)
        end if    
        call MPI_Bcast(n, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#else
        n = sqrt(n)
#endif        
    end function norm_complex_gridfun_tensorial_2D


   function norm_real_gridfun_tensorial_3D(this, u) result(n)
        class(grid_tensorial_3D) :: this
        real(kind=prec), intent(in), target :: u(:,:,:)
        real(kind=prec) :: n
#ifdef _OPENMP
        real(kind=prec), pointer :: v(:,:,:)
        real(kind=prec), pointer :: w(:)
        integer :: j
#endif           
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: n1
#endif        
#ifndef _OPENMP
        n = sum(spread(spread(this%weights_x, 2, this%n2max-this%n2min+1), 3, this%n3max-this%n3min+1) &
               *spread(spread(this%weights_y, 1, this%n1max-this%n1min+1), 3, this%n3max-this%n3min+1) &
               *spread(spread(this%weights_z, 1, this%n1max-this%n1min+1), 2, this%n2max-this%n2min+1) &
               *u**2)
#else
!$OMP PARALLEL DO PRIVATE(j, v, w) REDUCTION(+:n)
        do j=1,n_threads
            v => u(:,:,lbound(u,3)+this%jj(j-1):lbound(u,3)+this%jj(j)-1)        
            w => this%weights_z(lbound(this%weights_z,1)+this%jj(j-1):lbound(this%weights_z,1)+this%jj(j)-1)
            n = n + sum(spread(spread(this%weights_x, 2, this%n2max-this%n2min+1), 3, this%jj(j)-this%jj(j-1)) &
               *spread(spread(this%weights_y, 1, this%n1max-this%n1min+1), 3, this%jj(j)-this%jj(j-1) ) &
               *spread(spread(w, 1, this%n1max-this%n1min+1), 2, this%n2max-this%n2min+1) &
               *u**2)
        end do    
!$OMP END PARALLEL DO 
#endif
#ifdef _MPI_
        call MPI_Reduce(n, n1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            n = sqrt(n1)
        end if    
        call MPI_Bcast(n, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#else
        n = sqrt(n)
#endif        
    end function norm_real_gridfun_tensorial_3D


    function norm_complex_gridfun_tensorial_3D(this, u) result(n)
        class(grid_tensorial_3D) :: this
        complex(kind=prec), intent(in), target :: u(:,:,:)
        real(kind=prec) :: n
#ifdef _OPENMP
        complex(kind=prec), pointer :: v(:,:,:)
        real(kind=prec), pointer :: w(:)
        integer :: j
#endif               
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: n1
#endif        
#ifndef _OPENMP
        n = sum(spread(spread(this%weights_x, 2, this%n2max-this%n2min+1), 3, this%n3max-this%n3min+1) &
               *spread(spread(this%weights_y, 1, this%n1max-this%n1min+1), 3, this%n3max-this%n3min+1) &
               *spread(spread(this%weights_z, 1, this%n1max-this%n1min+1), 2, this%n2max-this%n2min+1) &
               *(real(u,kind=prec)**2 +aimag(u)**2))
#else
        n = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v, w) REDUCTION(+:n)
        do j=1,n_threads
            v => u(:,:,lbound(u,3)+this%jj(j-1):lbound(u,3)+this%jj(j)-1)        
            w => this%weights_z(lbound(this%weights_z,1)+this%jj(j-1):lbound(this%weights_z,1)+this%jj(j)-1)
            n = n + sum(spread(spread(this%weights_x, 2, this%n2max-this%n2min+1), 3, this%jj(j)-this%jj(j-1)) &
               *spread(spread(this%weights_y, 1, this%n1max-this%n1min+1), 3, this%jj(j)-this%jj(j-1) ) &
               *spread(spread(w, 1, this%n1max-this%n1min+1), 2, this%n2max-this%n2min+1) &
               *(real(v,kind=prec)**2 +aimag(v)**2))
        end do    
!$OMP END PARALLEL DO 
#endif
#ifdef _MPI_
        call MPI_Reduce(n, n1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            n = sqrt(n1)
        end if    
        call MPI_Bcast(n, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#else
        n = sqrt(n)
#endif        
    end function norm_complex_gridfun_tensorial_3D

!YYYYY

   subroutine set_real_gridfun_tensorial_1D(this, u, f)
        class(grid_tensorial_1D) :: this
        real(kind=prec), intent(inout) :: u(this%n1min:this%n1max)
        real(kind=prec), external :: f

        integer :: ix

!$OMP PARALLEL DO PRIVATE(ix) 
        do ix = this%n1min,this%n1max
            u(ix) = f(this%nodes_x(ix))
        end do
!$OMP END PARALLEL DO
    end subroutine set_real_gridfun_tensorial_1D

    subroutine set_complex_gridfun_tensorial_1D(this, u, f)
        class(grid_tensorial_1D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max)
        complex(kind=prec), external :: f

        integer :: ix

!$OMP PARALLEL DO PRIVATE(ix) 
        do ix = this%n1min,this%n1max
            u(ix) = f(this%nodes_x(ix))
        end do
!$OMP END PARALLEL DO
    end subroutine set_complex_gridfun_tensorial_1D

    subroutine rset_complex_gridfun_tensorial_1D(this, u, f)
        class(grid_tensorial_1D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max)
        real(kind=prec), external :: f

        integer :: ix

!$OMP PARALLEL DO PRIVATE(ix) 
        do ix = this%n1min,this%n1max
            u(ix) = cmplx(f(this%nodes_x(ix)), kind=prec)
        end do
!$OMP END PARALLEL DO
    end subroutine rset_complex_gridfun_tensorial_1D



    subroutine set_real_gridfun_tensorial_2D(this, u, f)
        class(grid_tensorial_2D) :: this
        real(kind=prec), intent(out) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
        real(kind=prec), external :: f

        integer :: ix, iy

!$OMP PARALLEL DO PRIVATE(ix, iy) 
        do iy = this%n2min, this%n2max
            do ix = this%n1min,this%n1max
                u(ix, iy) = f(this%nodes_x(ix), this%nodes_y(iy))
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine set_real_gridfun_tensorial_2D

    subroutine set_complex_gridfun_tensorial_2D(this, u, f)
        class(grid_tensorial_2D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
        complex(kind=prec), external :: f

        integer :: ix, iy

!$OMP PARALLEL DO PRIVATE(ix, iy) 
        do iy = this%n2min, this%n2max
            do ix = this%n1min,this%n1max
                u(ix, iy) = f(this%nodes_x(ix), this%nodes_y(iy))
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine set_complex_gridfun_tensorial_2D

    subroutine rset_complex_gridfun_tensorial_2D(this, u, f)
        class(grid_tensorial_2D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
        real(kind=prec), external :: f

        integer :: ix, iy

!$OMP PARALLEL DO PRIVATE(ix, iy)
        do iy = this%n2min, this%n2max
            do ix = this%n1min,this%n1max
                u(ix, iy) = cmplx(f(this%nodes_x(ix), this%nodes_y(iy)), kind=prec)
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine rset_complex_gridfun_tensorial_2D

    subroutine set_real_gridfun_tensorial_3D(this, u, f)
        class(grid_tensorial_3D) :: this
        real(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        real(kind=prec), external :: f

        integer :: ix, iy, iz

!$OMP PARALLEL DO PRIVATE(ix, iy, iz)
        do iz = this%n3min, this%n3max
            do iy = this%n2min, this%n2max
                do ix = this%n1min,this%n1max
                    u(ix, iy, iz) = f(this%nodes_x(ix), this%nodes_y(iy), this%nodes_z(iz))
                end do
            end do
        end do
!$OMP END PARALLEL DO
        
    end subroutine set_real_gridfun_tensorial_3D

    subroutine set_complex_gridfun_tensorial_3D(this, u, f)
        class(grid_tensorial_3D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        complex(kind=prec), external :: f

        integer :: ix, iy, iz

!$OMP PARALLEL DO PRIVATE(ix, iy, iz)
        do iz = this%n3min, this%n3max
            do iy = this%n2min, this%n2max
                do ix = this%n1min,this%n1max
                    u(ix, iy, iz) = f(this%nodes_x(ix), this%nodes_y(iy), this%nodes_z(iz))
                end do
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine set_complex_gridfun_tensorial_3D

    subroutine rset_complex_gridfun_tensorial_3D(this, u, f)
        class(grid_tensorial_3D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        real(kind=prec), external :: f

        integer :: ix, iy, iz

!$OMP PARALLEL DO PRIVATE(ix, iy, iz)
        do iz = this%n3min, this%n3max
            do iy = this%n2min, this%n2max
                do ix = this%n1min,this%n1max
                    u(ix, iy, iz) = cmplx(f(this%nodes_x(ix), this%nodes_y(iy), this%nodes_z(iz)), kind=prec)
                end do
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine rset_complex_gridfun_tensorial_3D


    subroutine compute_nodes_xy_polar_2D(this)
        class(grid_polar_2D) :: this
        integer :: itheta, ir
        real(kind=prec) :: r, theta
!$OMP PARALLEL DO PRIVATE(ir, itheta, r, theta)
        do itheta = this%n2min, this%n2max
            theta = itheta*this%dtheta 
            do ir = this%n1min,this%n1max
                r = this%nodes_r(ir)
                this%nodes_x(ir, itheta) = r*cos(theta)
                this%nodes_y(ir, itheta) = r*sin(theta)
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine compute_nodes_xy_polar_2D


    subroutine compute_nodes_xy_cylindrical_3D(this)
        class(grid_cylindrical_3D) :: this
        integer :: itheta, ir
        real(kind=prec) :: r, theta
!$OMP PARALLEL DO PRIVATE(ir, itheta, r, theta)
        do itheta = this%n3min, this%n3max
            theta = itheta*this%dtheta 
            do ir = this%n1min,this%n1max
                r = this%nodes_r(ir)
                this%nodes_x(ir, itheta) = r*cos(theta)
                this%nodes_y(ir, itheta) = r*sin(theta)
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine compute_nodes_xy_cylindrical_3D



    subroutine set_real_gridfun_polar_2D(this, u, f, polar_coordinates)
        class(grid_polar_2D) :: this
        real(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
        real(kind=prec), external :: f
        logical, intent(in), optional :: polar_coordinates
        
        integer :: itheta, ir
        real(kind=prec) :: r, theta, x, y

        if (present(polar_coordinates).and.polar_coordinates) then
!$OMP PARALLEL DO PRIVATE(ir, itheta, r, theta)
            do itheta = this%n2min, this%n2max
                theta = itheta*this%dtheta 
                do ir = this%n1min,this%n1max
                    r = this%nodes_r(ir)
                    u(ir, itheta) = cmplx(f(r, theta), kind=prec)
                end do
            end do
!$OMP END PARALLEL DO
        else
!$OMP PARALLEL DO PRIVATE(ir, itheta, r, theta, x, y)
            do itheta = this%n2min, this%n2max
        !       theta = itheta*this%dtheta 
                do ir = this%n1min,this%n1max
        !           r = this%nodes_r(ir)
        !           x = r*cos(theta)
        !           y = r*sin(theta)
        !           u(ir, itheta) = cmplx(f(x,y), kind=prec)
                    u(ir, itheta) = cmplx(f(this%nodes_x(ir,itheta),this%nodes_y(ir,itheta)), kind=prec)
                end do
            end do
!$OMP END PARALLEL DO
        end if
    end subroutine set_real_gridfun_polar_2D

    subroutine set_complex_gridfun_polar_2D(this, u, f, polar_coordinates)
        class(grid_polar_2D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
        complex(kind=prec), external :: f
        logical, intent(in), optional :: polar_coordinates
        
        integer :: itheta, ir
        real(kind=prec) :: r, theta, x, y

        if (present(polar_coordinates).and.polar_coordinates) then
!$OMP PARALLEL DO PRIVATE(ir, itheta, r, theta)
            do itheta = this%n2min, this%n2max
                theta = itheta*this%dtheta 
                do ir = this%n1min,this%n1max
                    r = this%nodes_r(ir)
                    u(ir, itheta) = cmplx(f(r, theta), kind=prec)
                end do
            end do
!$OMP END PARALLEL DO
        else
!$OMP PARALLEL DO PRIVATE(ir, itheta, r, theta, x, y)
            do itheta = this%n2min, this%n2max
        !       theta = itheta*this%dtheta 
                do ir = this%n1min,this%n1max
        !           r = this%nodes_r(ir)
        !           x = r*cos(theta)
        !           y = r*sin(theta)
        !           u(ir, itheta) = cmplx(f(x,y), kind=prec)
                    u(ir, itheta) = cmplx(f(this%nodes_x(ir,itheta),this%nodes_y(ir,itheta)), kind=prec)
                end do
            end do
!$OMP END PARALLEL DO
        end if
    end subroutine set_complex_gridfun_polar_2D

    subroutine rset_complex_gridfun_polar_2D(this, u, f, polar_coordinates)
        class(grid_polar_2D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
        real(kind=prec), external :: f
        logical, intent(in), optional :: polar_coordinates
        
        integer :: itheta, ir
        real(kind=prec) :: r, theta, x, y

        if (present(polar_coordinates).and.polar_coordinates) then
!$OMP PARALLEL DO PRIVATE(ir, itheta, r, theta)
            do itheta = this%n2min, this%n2max
                theta = itheta*this%dtheta 
                do ir = this%n1min,this%n1max
                    r = this%nodes_r(ir)
                    u(ir, itheta) = cmplx(f(r, theta), kind=prec)
                end do
            end do
!$OMP END PARALLEL DO
        else
!$OMP PARALLEL DO PRIVATE(ir, itheta, r, theta, x, y)
            do itheta = this%n2min, this%n2max
        !       theta = itheta*this%dtheta 
                do ir = this%n1min,this%n1max
        !           r = this%nodes_r(ir)
        !           x = r*cos(theta)
        !           y = r*sin(theta)
        !           u(ir, itheta) = cmplx(f(x,y), kind=prec)
                    u(ir, itheta) = cmplx(f(this%nodes_x(ir,itheta),this%nodes_y(ir,itheta)), kind=prec)
                end do
            end do
!$OMP END PARALLEL DO
        end if
    end subroutine rset_complex_gridfun_polar_2D
    

   subroutine set_real_gridfun_cylindrical_3D(this, u, f, cylindrical_coordinates)
        class(grid_cylindrical_3D) :: this
        real(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        real(kind=prec), external :: f
        logical, intent(in), optional :: cylindrical_coordinates
        
        integer :: itheta, ir, iz
        real(kind=prec) :: r, theta, x, y, z

        if (present(cylindrical_coordinates).and.cylindrical_coordinates) then
!$OMP PARALLEL DO PRIVATE(ir, itheta, iz, r, theta, z)
          do itheta = this%n3min, this%n3max
            theta = itheta*this%dtheta 
            do iz = this%n2min, this%n2max
                z = this%nodes_z(iz)
                do ir = this%n1min,this%n1max
                    r = this%nodes_r(ir)
                    theta = itheta*this%dtheta 
                    u(ir, iz, itheta) = cmplx(f(r, theta,z), kind=prec)
                end do
            end do
          end do
!$OMP END PARALLEL DO
        else
!$OMP PARALLEL DO PRIVATE(ir, itheta, iz, r, theta, x, y, z)
          do itheta = this%n3min, this%n3max
        !   theta = itheta*this%dtheta 
            do iz = this%n2min, this%n2max
                z = this%nodes_z(iz)
                do ir = this%n1min,this%n1max
         !          r = this%nodes_r(ir)
         !          x = r*cos(theta)
         !          y = r*sin(theta)
         !          u(ir, iz, itheta) = cmplx(f(x,y,z), kind=prec)
                    u(ir, iz, itheta) = cmplx(f(this%nodes_x(ir,itheta),this%nodes_y(ir,itheta),z), kind=prec)
                end do
            end do
          end do
!$OMP END PARALLEL DO
        end if
    end subroutine set_real_gridfun_cylindrical_3D

    subroutine set_complex_gridfun_cylindrical_3D(this, u, f, cylindrical_coordinates)
        class(grid_cylindrical_3D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        complex(kind=prec), external :: f
        logical, intent(in), optional :: cylindrical_coordinates
        
        integer :: itheta, ir, iz
        real(kind=prec) :: r, theta, x, y, z

        if (present(cylindrical_coordinates).and.cylindrical_coordinates) then
!$OMP PARALLEL DO PRIVATE(ir, itheta, iz, r, theta, z)
          do itheta = this%n3min, this%n3max
            theta = itheta*this%dtheta 
            do iz = this%n2min, this%n2max
                z = this%nodes_z(iz)
                do ir = this%n1min,this%n1max
                    r = this%nodes_r(ir)
                    u(ir, iz, itheta) = cmplx(f(r, theta,z), kind=prec)
                end do
            end do
          end do
!$OMP END PARALLEL DO
        else
!$OMP PARALLEL DO PRIVATE(ir, itheta, iz, r, theta, x, y, z)
          do itheta = this%n3min, this%n3max
        !   theta = itheta*this%dtheta 
            do iz = this%n2min, this%n2max
                z = this%nodes_z(iz)
                do ir = this%n1min,this%n1max
        !           r = this%nodes_r(ir)
        !           x = r*cos(theta)
        !           y = r*sin(theta)
        !           u(ir, iz, itheta) = cmplx(f(x,y,z), kind=prec)
                    u(ir, iz, itheta) = cmplx(f(this%nodes_x(ir,itheta),this%nodes_y(ir,itheta),z), kind=prec)
                end do
            end do
          end do
!$OMP END PARALLEL DO
        end if
    end subroutine set_complex_gridfun_cylindrical_3D

    subroutine rset_complex_gridfun_cylindrical_3D(this, u, f, cylindrical_coordinates)
        class(grid_cylindrical_3D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        real(kind=prec), external :: f
        logical, intent(in), optional :: cylindrical_coordinates
        
        integer :: itheta, ir, iz
        real(kind=prec) :: r, theta, x, y, z

        if (present(cylindrical_coordinates).and.cylindrical_coordinates) then
!$OMP PARALLEL DO PRIVATE(ir, itheta, iz, r, theta, z)
          do itheta = this%n3min, this%n3max
            theta = itheta*this%dtheta 
            do iz = this%n2min, this%n2max
                z = this%nodes_z(iz)
                do ir = this%n1min,this%n1max
                    r = this%nodes_r(ir)
                    u(ir, iz, itheta) = cmplx(f(r, theta, z), kind=prec)
                end do
            end do
          end do
!$OMP END PARALLEL DO
        else
!$OMP PARALLEL DO PRIVATE(ir, itheta, iz, r, theta, x, y, z)
          do itheta = this%n3min, this%n3max
        !   theta = itheta*this%dtheta 
            do iz = this%n2min, this%n2max
                z = this%nodes_z(iz)
                do ir = this%n1min,this%n1max
        !           r = this%nodes_r(ir)
        !           x = r*cos(theta)
        !           y = r*sin(theta)
        !           u(ir, iz, itheta) = cmplx(f(x,y,z), kind=prec)
                    u(ir, iz, itheta) = cmplx(f(this%nodes_x(ir,itheta),this%nodes_y(ir,itheta),z), kind=prec)
                end do
            end do
          end do
!$OMP END PARALLEL DO
        end if
    end subroutine rset_complex_gridfun_cylindrical_3D

!###
!YYYYY

   subroutine set_t_real_gridfun_tensorial_1D(this, u, f, t)
        class(grid_tensorial_1D) :: this
        real(kind=prec), intent(inout) :: u(this%n1min:this%n1max)
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t

        integer :: ix

!$OMP PARALLEL DO PRIVATE(ix) 
        do ix = this%n1min,this%n1max
            u(ix) = f(this%nodes_x(ix), t)
        end do
!$OMP END PARALLEL DO
    end subroutine set_t_real_gridfun_tensorial_1D

    subroutine set_t_complex_gridfun_tensorial_1D(this, u, f, t)
        class(grid_tensorial_1D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max)
        complex(kind=prec), external :: f
        real(kind=prec), intent(in) :: t

        integer :: ix

!$OMP PARALLEL DO PRIVATE(ix) 
        do ix = this%n1min,this%n1max
            u(ix) = f(this%nodes_x(ix), t)
        end do
!$OMP END PARALLEL DO
    end subroutine set_t_complex_gridfun_tensorial_1D

    subroutine rset_t_complex_gridfun_tensorial_1D(this, u, f, t)
        class(grid_tensorial_1D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max)
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t

        integer :: ix

!$OMP PARALLEL DO PRIVATE(ix) 
        do ix = this%n1min,this%n1max
            u(ix) = cmplx(f(this%nodes_x(ix), t), kind=prec)
        end do
!$OMP END PARALLEL DO
    end subroutine rset_t_complex_gridfun_tensorial_1D



    subroutine set_t_real_gridfun_tensorial_2D(this, u, f, t)
        class(grid_tensorial_2D) :: this
        real(kind=prec), intent(out) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t

        integer :: ix, iy

!$OMP PARALLEL DO PRIVATE(ix, iy) 
        do iy = this%n2min, this%n2max
            do ix = this%n1min,this%n1max
                u(ix, iy) = f(this%nodes_x(ix), this%nodes_y(iy), t)
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine set_t_real_gridfun_tensorial_2D

    subroutine set_t_complex_gridfun_tensorial_2D(this, u, f, t)
        class(grid_tensorial_2D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
        complex(kind=prec), external :: f
        real(kind=prec), intent(in) :: t

        integer :: ix, iy

!$OMP PARALLEL DO PRIVATE(ix, iy) 
        do iy = this%n2min, this%n2max
            do ix = this%n1min,this%n1max
                u(ix, iy) = f(this%nodes_x(ix), this%nodes_y(iy), t)
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine set_t_complex_gridfun_tensorial_2D

    subroutine rset_t_complex_gridfun_tensorial_2D(this, u, f, t)
        class(grid_tensorial_2D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t

        integer :: ix, iy

!$OMP PARALLEL DO PRIVATE(ix, iy)
        do iy = this%n2min, this%n2max
            do ix = this%n1min,this%n1max
                u(ix, iy) = cmplx(f(this%nodes_x(ix), this%nodes_y(iy), t), kind=prec)
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine rset_t_complex_gridfun_tensorial_2D

    subroutine set_t_real_gridfun_tensorial_3D(this, u, f, t)
        class(grid_tensorial_3D) :: this
        real(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t

        integer :: ix, iy, iz

!$OMP PARALLEL DO PRIVATE(ix, iy, iz)
        do iz = this%n3min, this%n3max
            do iy = this%n2min, this%n2max
                do ix = this%n1min,this%n1max
                    u(ix, iy, iz) = f(this%nodes_x(ix), this%nodes_y(iy), this%nodes_z(iz), t)
                end do
            end do
        end do
!$OMP END PARALLEL DO
        
    end subroutine set_t_real_gridfun_tensorial_3D

    subroutine set_t_complex_gridfun_tensorial_3D(this, u, f, t)
        class(grid_tensorial_3D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        complex(kind=prec), external :: f
        real(kind=prec), intent(in) :: t

        integer :: ix, iy, iz

!$OMP PARALLEL DO PRIVATE(ix, iy, iz)
        do iz = this%n3min, this%n3max
            do iy = this%n2min, this%n2max
                do ix = this%n1min,this%n1max
                    u(ix, iy, iz) = f(this%nodes_z(ix), this%nodes_y(iy), this%nodes_z(iz), t)
                end do
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine set_t_complex_gridfun_tensorial_3D

    subroutine rset_t_complex_gridfun_tensorial_3D(this, u, f, t)
        class(grid_tensorial_3D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t

        integer :: ix, iy, iz

!$OMP PARALLEL DO PRIVATE(ix, iy, iz)
        do iz = this%n3min, this%n3max
            do iy = this%n2min, this%n2max
                do ix = this%n1min,this%n1max
                    u(ix, iy, iz) = cmplx(f(this%nodes_z(ix), this%nodes_y(iy), this%nodes_z(iz), t), kind=prec)
                end do
            end do
        end do
!$OMP END PARALLEL DO
    end subroutine rset_t_complex_gridfun_tensorial_3D


    subroutine set_t_real_gridfun_polar_2D(this, u, f, t, polar_coordinates)
        class(grid_polar_2D) :: this
        real(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t
        logical, intent(in), optional :: polar_coordinates
        
        integer :: itheta, ir
        real(kind=prec) :: r, theta, x, y

        if (present(polar_coordinates).and.polar_coordinates) then
!$OMP PARALLEL DO PRIVATE(ir, itheta, r, theta)
            do itheta = this%n2min, this%n2max
                theta = itheta*this%dtheta 
                do ir = this%n1min,this%n1max
                    r = this%nodes_r(ir)
                    u(ir, itheta) = cmplx(f(r, theta, t), kind=prec)
                end do
            end do
!$OMP END PARALLEL DO
        else
!$OMP PARALLEL DO PRIVATE(ir, itheta, r, theta, x, y)
            do itheta = this%n2min, this%n2max
                theta = itheta*this%dtheta 
                do ir = this%n1min,this%n1max
                    r = this%nodes_r(ir)
                    x = r*cos(theta)
                    y = r*sin(theta)
                    u(ir, itheta) = cmplx(f(x,y, t), kind=prec)
                end do
            end do
!$OMP END PARALLEL DO
        end if
    end subroutine set_t_real_gridfun_polar_2D

    subroutine set_t_complex_gridfun_polar_2D(this, u, f, t, polar_coordinates)
        class(grid_polar_2D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
        complex(kind=prec), external :: f
        real(kind=prec), intent(in) :: t
        logical, intent(in), optional :: polar_coordinates
        
        integer :: itheta, ir
        real(kind=prec) :: r, theta, x, y

        if (present(polar_coordinates).and.polar_coordinates) then
!$OMP PARALLEL DO PRIVATE(ir, itheta, r, theta)
            do itheta = this%n2min, this%n2max
                theta = itheta*this%dtheta 
                do ir = this%n1min,this%n1max
                    r = this%nodes_r(ir)
                    u(ir, itheta) = cmplx(f(r, theta, t), kind=prec)
                end do
            end do
!$OMP END PARALLEL DO
        else
!$OMP PARALLEL DO PRIVATE(ir, itheta, r, theta, x, y)
            do itheta = this%n2min, this%n2max
                theta = itheta*this%dtheta 
                do ir = this%n1min,this%n1max
                    r = this%nodes_r(ir)
                    x = r*cos(theta)
                    y = r*sin(theta)
                    u(ir, itheta) = cmplx(f(x,y, t), kind=prec)
                end do
            end do
!$OMP END PARALLEL DO
        end if
    end subroutine set_t_complex_gridfun_polar_2D

    subroutine rset_t_complex_gridfun_polar_2D(this, u, f, t, polar_coordinates)
        class(grid_polar_2D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max)
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t
        logical, intent(in), optional :: polar_coordinates
        
        integer :: itheta, ir
        real(kind=prec) :: r, theta, x, y

        if (present(polar_coordinates).and.polar_coordinates) then
!$OMP PARALLEL DO PRIVATE(ir, itheta, r, theta)
            do itheta = this%n2min, this%n2max
                theta = itheta*this%dtheta 
                do ir = this%n1min,this%n1max
                    r = this%nodes_r(ir)
                    u(ir, itheta) = cmplx(f(r, theta, t), kind=prec)
                end do
            end do
!$OMP END PARALLEL DO
        else
!$OMP PARALLEL DO PRIVATE(ir, itheta, r, theta, x, y)
            do itheta = this%n2min, this%n2max
                theta = itheta*this%dtheta 
                do ir = this%n1min,this%n1max
                    r = this%nodes_r(ir)
                    x = r*cos(theta)
                    y = r*sin(theta)
                    u(ir, itheta) = cmplx(f(x,y, t), kind=prec)
                end do
            end do
!$OMP END PARALLEL DO
        end if
    end subroutine rset_t_complex_gridfun_polar_2D
    

   subroutine set_t_real_gridfun_cylindrical_3D(this, u, f, t, cylindrical_coordinates)
        class(grid_cylindrical_3D) :: this
        real(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t
        logical, intent(in), optional :: cylindrical_coordinates
        
        integer :: itheta, ir, iz
        real(kind=prec) :: r, theta, x, y, z

        if (present(cylindrical_coordinates).and.cylindrical_coordinates) then
!$OMP PARALLEL DO PRIVATE(ir, itheta, iz, r, theta, z)
          do itheta = this%n3min, this%n3max
            theta = itheta*this%dtheta 
            do iz = this%n2min, this%n2max
                z = this%nodes_z(iz)
                do ir = this%n1min,this%n1max
                    r = this%nodes_r(ir)
                    theta = itheta*this%dtheta 
                    u(ir, iz, itheta) = cmplx(f(r, theta,z, t), kind=prec)
                end do
            end do
          end do
!$OMP END PARALLEL DO
        else
!$OMP PARALLEL DO PRIVATE(ir, itheta, iz, r, theta, x, y, z)
          do itheta = this%n3min, this%n3max
            theta = itheta*this%dtheta 
            do iz = this%n2min, this%n2max
                z = this%nodes_z(iz)
                do ir = this%n1min,this%n1max
                    r = this%nodes_r(ir)
                    x = r*cos(theta)
                    y = r*sin(theta)
                    u(ir, iz, itheta) = cmplx(f(x,y,z, t), kind=prec)
                end do
            end do
          end do
!$OMP END PARALLEL DO
        end if
    end subroutine set_t_real_gridfun_cylindrical_3D

    subroutine set_t_complex_gridfun_cylindrical_3D(this, u, f, t, cylindrical_coordinates)
        class(grid_cylindrical_3D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        complex(kind=prec), external :: f
        real(kind=prec), intent(in) :: t
        logical, intent(in), optional :: cylindrical_coordinates
        
        integer :: itheta, ir, iz
        real(kind=prec) :: r, theta, x, y, z

        if (present(cylindrical_coordinates).and.cylindrical_coordinates) then
!$OMP PARALLEL DO PRIVATE(ir, itheta, iz, r, theta, z)
          do itheta = this%n3min, this%n3max
            theta = itheta*this%dtheta 
            do iz = this%n2min, this%n2max
                z = this%nodes_z(iz)
                do ir = this%n1min,this%n1max
                    r = this%nodes_r(ir)
                    u(ir, iz, itheta) = cmplx(f(r, theta,z, t), kind=prec)
                end do
            end do
          end do
!$OMP END PARALLEL DO
        else
!$OMP PARALLEL DO PRIVATE(ir, itheta, iz, r, theta, x, y, z)
          do itheta = this%n3min, this%n3max
            theta = itheta*this%dtheta 
            do iz = this%n2min, this%n2max
                z = this%nodes_z(iz)
                do ir = this%n1min,this%n1max
                    r = this%nodes_r(ir)
                    x = r*cos(theta)
                    y = r*sin(theta)
                    u(ir, iz, itheta) = cmplx(f(x,y,z, t), kind=prec)
                end do
            end do
          end do
!$OMP END PARALLEL DO
        end if
    end subroutine set_t_complex_gridfun_cylindrical_3D

    subroutine rset_t_complex_gridfun_cylindrical_3D(this, u, f, t, cylindrical_coordinates)
        class(grid_cylindrical_3D) :: this
        complex(kind=prec), intent(inout) :: u(this%n1min:this%n1max,this%n2min:this%n2max,this%n3min:this%n3max)
        real(kind=prec), external :: f
        real(kind=prec), intent(in) :: t
        logical, intent(in), optional :: cylindrical_coordinates
        
        integer :: itheta, ir, iz
        real(kind=prec) :: r, theta, x, y, z

        if (present(cylindrical_coordinates).and.cylindrical_coordinates) then
!$OMP PARALLEL DO PRIVATE(ir, itheta, iz, r, theta, z)
          do itheta = this%n3min, this%n3max
            theta = itheta*this%dtheta 
            do iz = this%n2min, this%n2max
                z = this%nodes_z(iz)
                do ir = this%n1min,this%n1max
                    r = this%nodes_r(ir)
                    u(ir, iz, itheta) = cmplx(f(r, theta, z, t), kind=prec)
                end do
            end do
          end do
!$OMP END PARALLEL DO
        else
!$OMP PARALLEL DO PRIVATE(ir, itheta, iz, r, theta, x, y, z)
          do itheta = this%n3min, this%n3max
            theta = itheta*this%dtheta 
            do iz = this%n2min, this%n2max
                z = this%nodes_z(iz)
                do ir = this%n1min,this%n1max
                    r = this%nodes_r(ir)
                    x = r*cos(theta)
                    y = r*sin(theta)
                    u(ir, iz, itheta) = cmplx(f(x,y,z, t), kind=prec)
                end do
            end do
          end do
!$OMP END PARALLEL DO
        end if
    end subroutine rset_t_complex_gridfun_cylindrical_3D

!###


   function norm_real_gridfun_polar_2D(this, u) result(n)
        class(grid_polar_2D) :: this
        real(kind=prec), intent(in), target :: u(:,:)
        real(kind=prec) :: n
#ifdef _OPENMP
        real(kind=prec), pointer :: v(:,:)
        integer :: j
#endif                  
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: n1
#endif        
#ifndef _OPENMP
        n = sum( spread(this%weights_r, 2, this%n2max-this%n2min+1) *u**2 )*this%dtheta
#else
        n = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v) REDUCTION(+:n)
        do j=1,n_threads
            v => u(:,lbound(u,2)+this%jj(j-1):lbound(u,2)+this%jj(j)-1)        
            n = n + sum(spread(this%weights_r, 2, this%jj(j)-this%jj(j-1))*v**2)*this%dtheta
        end do    
!$OMP END PARALLEL DO 
#endif 
#ifdef _MPI_
        call MPI_Reduce(n, n1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            n = sqrt(n1)
        end if    
        call MPI_Bcast(n, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#else
        n = sqrt(n)
#endif        
    end function norm_real_gridfun_polar_2D


   function inner_product_real_gridfun_polar_2D(this, u1, u2) result(n)
        class(grid_polar_2D) :: this
        real(kind=prec), intent(in), target :: u1(:,:)
        real(kind=prec), intent(in), target :: u2(:,:)
        real(kind=prec) :: n
#ifdef _OPENMP
        real(kind=prec), pointer :: v1(:,:)
        real(kind=prec), pointer :: v2(:,:)
        integer :: j
#endif                  
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: n1
#endif        
#ifndef _OPENMP
        n = sum( spread(this%weights_r, 2, this%n2max-this%n2min+1) *u1*u2 )*this%dtheta
#else
        n = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v1, v2) REDUCTION(+:n)
        do j=1,n_threads
            v1 => u1(:,lbound(u1,2)+this%jj(j-1):lbound(u1,2)+this%jj(j)-1)        
            v2 => u2(:,lbound(u2,2)+this%jj(j-1):lbound(u2,2)+this%jj(j)-1)        
            n = n + sum(spread(this%weights_r, 2, this%jj(j)-this%jj(j-1))*v1*v2)*this%dtheta
        end do    
!$OMP END PARALLEL DO 
#endif 
#ifdef _MPI_
        call MPI_Reduce(n, n1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            n = n1
        end if    
        call MPI_Bcast(n, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#else
        n = n
#endif        
    end function inner_product_real_gridfun_polar_2D




    function norm_complex_gridfun_polar_2D(this, u) result(n)
        class(grid_polar_2D) :: this
        complex(kind=prec), intent(in), target :: u(:,:)
        real(kind=prec) :: n
#ifdef _OPENMP
        complex(kind=prec), pointer :: v(:,:)
        integer :: j
#endif                    
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: n1
#endif  
#ifndef _OPENMP
        n = sum(spread(this%weights_r, 2, this%n2max-this%n2min+1) &
               *(real(u,kind=prec)**2 +aimag(u)**2)) * this%dtheta
#else
        n = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v) REDUCTION(+:n)
        do j=1,n_threads
            v => u(:,lbound(u,2)+this%jj(j-1):lbound(u,2)+this%jj(j)-1)        
            n = n + sum(spread(this%weights_r, 2, this%jj(j)-this%jj(j-1) ) &
                       *(real(v,kind=prec)**2 +aimag(v)**2)) * this%dtheta
        end do    
!$OMP END PARALLEL DO 
#endif
#ifdef _MPI_
        call MPI_Reduce(n, n1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            n = sqrt(n1)
        end if    
        call MPI_Bcast(n, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#else
        n = sqrt(n)
#endif        
    end function norm_complex_gridfun_polar_2D

    function inner_product_complex_gridfun_polar_2D(this, u1, u2) result(n)
        class(grid_polar_2D) :: this
        complex(kind=prec), intent(in), target :: u1(:,:)
        complex(kind=prec), intent(in), target :: u2(:,:)
        complex(kind=prec) :: n
#ifdef _OPENMP
        complex(kind=prec), pointer :: v1(:,:)
        complex(kind=prec), pointer :: v2(:,:)
        integer :: j
#endif                    
#ifdef _MPI_
        integer :: ierr
        complex(kind=prec) :: n1
#endif  
#ifndef _OPENMP
        n = sum(spread(this%weights_r, 2, this%n2max-this%n2min+1) &
               *(conjg(u1)*u2)) * this%dtheta
#else
        n = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v1, v2) REDUCTION(+:n)
        do j=1,n_threads
            v1 => u1(:,lbound(u1,2)+this%jj(j-1):lbound(u1,2)+this%jj(j)-1)        
            v2 => u2(:,lbound(u2,2)+this%jj(j-1):lbound(u2,2)+this%jj(j)-1)        
            n = n + sum(spread(this%weights_r, 2, this%jj(j)-this%jj(j-1) ) &
                       *(conjg(v1)*v2)) * this%dtheta
        end do    
!$OMP END PARALLEL DO 
#endif
#ifdef _MPI_
        call MPI_Reduce(n, n1, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            n = n1
        end if    
        call MPI_Bcast(n, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr) 
#else
        n = n
#endif        
    end function inner_product_complex_gridfun_polar_2D


  
   function norm_real_gridfun_cylindrical_3D(this, u) result(n)
        class(grid_cylindrical_3D) :: this
        real(kind=prec), intent(in), target :: u(:,:,:)
        real(kind=prec) :: n
#ifdef _OPENMP
        real(kind=prec), pointer :: v(:,:,:)
        integer :: j
#endif           
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: n1
#endif        
#ifndef _OPENMP
        n = sum(spread(spread(this%weights_r, 2, this%n2max-this%n2min+1), 3, this%n3max-this%n3min+1) &
               *spread(spread(this%weights_z, 1, this%n1max-this%n1min+1), 3, this%n3max-this%n3min+1) &
               *u**2) * this%dtheta
#else
!$OMP PARALLEL DO PRIVATE(j, v) REDUCTION(+:n)
        do j=1,n_threads
            v => u(:,:,lbound(u,3)+this%jj(j-1):lbound(u,3)+this%jj(j)-1)        
            n = n + sum(spread(spread(this%weights_r, 2, this%n2max-this%n2min+1), 3, this%jj(j)-this%jj(j-1)) &
               *spread(spread(this%weights_z, 1, this%n1max-this%n1min+1), 3, this%jj(j)-this%jj(j-1) ) &
               *v**2) + this%dtheta
        end do    
!$OMP END PARALLEL DO 
#endif
#ifdef _MPI_
        call MPI_Reduce(n, n1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            n = sqrt(n1)
        end if    
        call MPI_Bcast(n, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#else
        n = sqrt(n)
#endif        
    end function norm_real_gridfun_cylindrical_3D

   function inner_product_real_gridfun_cylindrical_3D(this, u1, u2) result(n)
        class(grid_cylindrical_3D) :: this
        real(kind=prec), intent(in), target :: u1(:,:,:)
        real(kind=prec), intent(in), target :: u2(:,:,:)
        real(kind=prec) :: n
#ifdef _OPENMP
        real(kind=prec), pointer :: v1(:,:,:)
        real(kind=prec), pointer :: v2(:,:,:)
        integer :: j
#endif           
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: n1
#endif        
#ifndef _OPENMP
        n = sum(spread(spread(this%weights_r, 2, this%n2max-this%n2min+1), 3, this%n3max-this%n3min+1) &
               *spread(spread(this%weights_z, 1, this%n1max-this%n1min+1), 3, this%n3max-this%n3min+1) &
               *u1*u2) * this%dtheta
#else
!$OMP PARALLEL DO PRIVATE(j, v1, v2) REDUCTION(+:n)
        do j=1,n_threads
            v1 => u1(:,:,lbound(u1,3)+this%jj(j-1):lbound(u1,3)+this%jj(j)-1)        
            v2 => u2(:,:,lbound(u2,3)+this%jj(j-1):lbound(u2,3)+this%jj(j)-1)        
            n = n + sum(spread(spread(this%weights_r, 2, this%n2max-this%n2min+1), 3, this%jj(j)-this%jj(j-1)) &
               *spread(spread(this%weights_z, 1, this%n1max-this%n1min+1), 3, this%jj(j)-this%jj(j-1) ) &
               *v1*v2) + this%dtheta
        end do    
!$OMP END PARALLEL DO 
#endif
#ifdef _MPI_
        call MPI_Reduce(n, n1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            n = n1
        end if    
        call MPI_Bcast(n, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#else
        n = n
#endif        
    end function inner_product_real_gridfun_cylindrical_3D



    function norm_complex_gridfun_cylindrical_3D(this, u) result(n)
        class(grid_cylindrical_3D) :: this
        complex(kind=prec), intent(in), target :: u(:,:,:)
        real(kind=prec) :: n
#ifdef _OPENMP
        complex(kind=prec), pointer :: v(:,:,:)
        integer :: j
#endif               
#ifdef _MPI_
        integer :: ierr
        real(kind=prec) :: n1
#endif        
#ifndef _OPENMP
        n = sum(spread(spread(this%weights_r, 2, this%n2max-this%n2min+1), 3, this%n3max-this%n3min+1) &
               *spread(spread(this%weights_z, 1, this%n1max-this%n1min+1), 3, this%n3max-this%n3min+1) &
               *(real(u,kind=prec)**2 +aimag(u)**2)) * this%dtheta
#else
        n = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v) REDUCTION(+:n)
        do j=1,n_threads
            v => u(:,:,lbound(u,3)+this%jj(j-1):lbound(u,3)+this%jj(j)-1)        
            n = n + sum(spread(spread(this%weights_r, 2, this%n2max-this%n2min+1), 3, this%jj(j)-this%jj(j-1)) &
               *spread(spread(this%weights_z, 1, this%n1max-this%n1min+1), 3, this%jj(j)-this%jj(j-1) ) &
               *(real(v,kind=prec)**2 +aimag(v)**2)) * this%dtheta
        end do    
!$OMP END PARALLEL DO 
#endif
#ifdef _MPI_
        call MPI_Reduce(n, n1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            n = sqrt(n1)
        end if    
        call MPI_Bcast(n, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#else
        n = sqrt(n)
#endif        
    end function norm_complex_gridfun_cylindrical_3D

  
    function inner_product_complex_gridfun_cylindrical_3D(this, u1, u2) result(n)
        class(grid_cylindrical_3D) :: this
        complex(kind=prec), intent(in), target :: u1(:,:,:)
        complex(kind=prec), intent(in), target :: u2(:,:,:)
        complex(kind=prec) :: n
#ifdef _OPENMP
        complex(kind=prec), pointer :: v1(:,:,:)
        complex(kind=prec), pointer :: v2(:,:,:)
        integer :: j
#endif               
#ifdef _MPI_
        integer :: ierr
        complex(kind=prec) :: n1
#endif        
#ifndef _OPENMP
        n = sum(spread(spread(this%weights_r, 2, this%n2max-this%n2min+1), 3, this%n3max-this%n3min+1) &
               *spread(spread(this%weights_z, 1, this%n1max-this%n1min+1), 3, this%n3max-this%n3min+1) &
               *(conjg(u1)*u2)) * this%dtheta
#else
        n = 0.0_prec
!$OMP PARALLEL DO PRIVATE(j, v1, v2) REDUCTION(+:n)
        do j=1,n_threads
            v1 => u1(:,:,lbound(u1,3)+this%jj(j-1):lbound(u1,3)+this%jj(j)-1)        
            v2 => u2(:,:,lbound(u2,3)+this%jj(j-1):lbound(u2,3)+this%jj(j)-1)        
            n = n + sum(spread(spread(this%weights_r, 2, this%n2max-this%n2min+1), 3, this%jj(j)-this%jj(j-1)) &
               *spread(spread(this%weights_z, 1, this%n1max-this%n1min+1), 3, this%jj(j)-this%jj(j-1) ) &
               *(conjg(v1)*v2)) * this%dtheta
        end do    
!$OMP END PARALLEL DO 
#endif
#ifdef _MPI_
        call MPI_Reduce(n, n1, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (this_proc==0) then
            n = n1
        end if    
        call MPI_Bcast(n, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr) 
#else
        n = n
#endif        
    end function inner_product_complex_gridfun_cylindrical_3D


#ifdef _QUADPRECISION_
end module tssmq_grid
#else
end module tssm_grid
#endif
