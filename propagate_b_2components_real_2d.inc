!-----------------------------------------------------------------------
! propagate_b_2components_real_2d.inc
!-----------------------------------------------------------------------

    elemental subroutine Explicit_Euler_1st_order(dt, u, v, u1, v1) 
        real(kind=prec), intent(in) :: dt
        real(kind=prec), intent(in) :: u, v
        real(kind=prec), intent(out) :: u1, v1
        real(kind=prec) :: k1
        real(kind=prec) :: l1
        k1 = F1(u,v)
        l1 = F2(u,v)
        u1 = u + dt*k1
        v1 = v + dt*l1
    end  subroutine Explicit_Euler_1st_order

    elemental subroutine Explicit_midpoint_2nd_order(dt, u, v, u1, v1) 
        real(kind=prec), intent(in) :: dt
        real(kind=prec), intent(in) :: u, v
        real(kind=prec), intent(out) :: u1, v1
        real(kind=prec) :: k1, k2
        real(kind=prec) :: l1, l2
        k1 = F1(u,v)
        l1 = F2(u,v)
        k2 = F1(u+dt*0.5_prec*k1, v+dt*0.5_prec*l1)
        l2 = F2(u+dt*0.5_prec*k1, v+dt*0.5_prec*l1)
        u1 = u + dt*k2
        v1 = v + dt*l2
    end  subroutine Explicit_midpoint_2nd_order

    elemental subroutine Heun_3rd_order(dt, u, v, u1, v1) 
        real(kind=prec), intent(in) :: dt
        real(kind=prec), intent(in) :: u, v
        real(kind=prec), intent(out) :: u1, v1
        real(kind=prec) :: k1, k2, k3 
        real(kind=prec) :: l1, l2, l3 
        k1 = F1(u,v)
        l1 = F2(u,v)
        k2 = F1(u+dt/3.0_prec*k1, v+dt/3.0_prec*l1)
        l2 = F2(u+dt/3.0_prec*k1, v+dt/3.0_prec*l1)
        k3 = F1(u+dt*2.0_prec/3.0_prec*k2, v+dt*2.0_prec/3.0_prec*l2)
        l3 = F2(u+dt*2.0_prec/3.0_prec*k2, v+dt*2.0_prec/3.0_prec*l2)
        u1 = u+dt*(0.25_prec*k1+0.75_prec*k3)
        v1 = v+dt*(0.25_prec*l1+0.75_prec*l3)
    end subroutine Heun_3rd_order

    elemental subroutine Runge_Kutta_4th_order(dt, u, v, u1, v1) 
        real(kind=prec), intent(in) :: dt
        real(kind=prec), intent(in) :: u, v
        real(kind=prec), intent(out) :: u1, v1
        real(kind=prec) :: k1, k2, k3, k4 
        real(kind=prec) :: l1, l2, l3, l4 
        k1 = F1(u, v)
        l1 = F2(u, v)
        k2 = F1(u+dt*0.5_prec*k1, v+dt*0.5_prec*l1)
        l2 = F2(u+dt*0.5_prec*k1, v+dt*0.5_prec*l1)
        k3 = F1(u+dt*0.5_prec*k2, v+dt*0.5_prec*l2)
        l3 = F2(u+dt*0.5_prec*k2, v+dt*0.5_prec*l2)
        k4 = F1(u+dt*k3, v+dt*l3)
        l4 = F2(u+dt*k3, v+dt*l3)
        u1 = u+dt*(1.0_prec/6.0_prec*k1+2.0_prec/6.0_prec*k2 &
                  +2.0_prec/6.0_prec*k3+1.0_prec/6.0_prec*k4)        
        v1 = v+dt*(1.0_prec/6.0_prec*l1+2.0_prec/6.0_prec*l2 &
                  +2.0_prec/6.0_prec*l3+1.0_prec/6.0_prec*l4)        
    end subroutine Runge_Kutta_4th_order

    subroutine INCLUDED_PROPAGATE_B(this, dt)
       class(wf_multicomponent_fourier_real_2D), intent(inout) :: this
       real(kind=prec), intent(in) :: dt

       real(kind=prec), dimension(:,:), pointer :: u1, u2

       call this%to_real_space  

       u1 => this%C(1)%p%u
       u2 => this%C(2)%p%u

       select case(this%m%method_for_B)
       case(1) 
          call Explicit_Euler_1st_order(dt, u1, u2, u1, u2)
       case(2) 
         call Explicit_midpoint_2nd_order(dt, u1, u2, u1, u2)
       case(3)                                         
          call Heun_3rd_order(dt, u1, u2, u1, u2)
       case(4)
          call Runge_Kutta_4th_order(dt, u1, u2, u1, u2)
       case default
          stop "E: method for B not implemented"
       end select   

   end subroutine INCLUDED_PROPAGATE_B

    !-----------------------------------------------------------------------

