module tssm_gray_scott    
    use tssm_fourier
    implicit none

    !parameter
    real(kind=prec), parameter :: D1 = 8E-5_prec
    real(kind=prec), parameter :: D2 = 4E-5_prec
    real(kind=prec), parameter :: gamma = 0.024_prec
    real(kind=prec), parameter :: kappa = 0.06_prec

    type, extends(multicomponent_fourier_2D) :: gray_scott
    end type gray_scott

    interface gray_scott
        module procedure :: new_gray_scott
    end interface gray_scott

    type, extends(wf_multicomponent_fourier_2D) :: wf_gray_scott
    contains 
        procedure :: propagate_B => propagate_B_wf_gray_scott
    end type wf_gray_scott

    interface wf_gray_scott
        module procedure new_wf_gray_scott
    end interface wf_gray_scott

contains
    function new_gray_scott(nx, xmin, xmax, ny, ymin, ymax, boundary_conditions) result(this)
        type(gray_scott) :: this
        real(kind=prec), intent(in) :: xmin 
        real(kind=prec), intent(in) :: xmax
        integer, intent(in) :: nx
        real(kind=prec), intent(in) :: ymin 
        real(kind=prec), intent(in) :: ymax
        integer, intent(in) :: ny
        integer, intent(in), optional :: boundary_conditions

        this%multicomponent_fourier_2D = multicomponent_fourier_2D(  &
                           (/ cmplx(D1, 0.0_prec, kind=prec),  cmplx(D2, 0.0_prec, kind=prec) /), &
                           nx, xmin, xmax, ny, ymin, ymax, boundary_conditions)
    end function new_gray_scott

    function new_wf_gray_scott(m) result(this)
        type(wf_gray_scott) :: this
        class(gray_scott), target, intent(inout) :: m
        
        this%wf_multicomponent_fourier_2D = wf_multicomponent_fourier_2D(m)
    end function new_wf_gray_scott

    !-----------------------------------------------------------------------

    elemental function F1(u,v)
        complex(kind=prec), intent(in) :: u
        complex(kind=prec), intent(in) :: v
        complex(kind=prec) :: F1
        F1 = -u*v**2 + gamma*(1.0_prec-u)
    end function F1

    elemental function F2(u,v)
        complex(kind=prec), intent(in) :: u
        complex(kind=prec), intent(in) :: v
        complex(kind=prec) :: F2
        F2 = u*v**2 - (gamma+kappa)*v
    end function F2

    include "propagate_b_2components_complex_2d.inc"

    subroutine propagate_B_wf_gray_scott(this, dt)
       class(wf_gray_scott), intent(inout) :: this
       complex(kind=prec), intent(in) :: dt
       call INCLUDED_PROPAGATE_B(this, dt)
    end subroutine propagate_B_wf_gray_scott


    !-----------------------------------------------------------------------
 

    function initial_value(x, y, which_component) result(u)
        real(kind=prec) :: u
        real(kind=prec), intent(in) :: x, y
        integer, intent(in) :: which_component

        if((1.0_prec<=x).and.(x<=1.5_prec).and. &
           (1.0_prec<=y).and.(y<=1.5_prec)) then
           u = 0.25_prec*(sin(4.0_prec*pi*x)*sin(4.0_prec*pi*y))**2
        else   
           u = 0.0_prec
        end if   
        if (which_component==1) then
           u = 1.0_prec - u
        end if   
    end function initial_value

    subroutine solution_out(psi, t, flag, abort_flag)
        class(wave_function), intent(inout) :: psi
        real(kind=prec), intent(in) :: t 
        character, intent(in) :: flag
        logical, intent(out) :: abort_flag
        integer, save :: step
        real(kind=prec) :: norm_U, norm_V

        select type (psi); class is (wf_gray_scott)

        if (flag=="i" .or. flag=="I") then ! init =======================
            step = 0 
            !if (this_proc==0) then
                write (*,*), "                       t               ||U||               " &
                             // "||V||"
                write (*,*), "-----------------------------------------------------------" &
                             //"------"
            !end if   
        end if   

        if (flag=="e" .or. flag=="E") then ! execute =====
            step = step + 1
               norm_U = psi%C(1)%norm2()
               norm_V = psi%C(2)%norm2()
             !  if (this_proc==0) then
                   write (*, '(I5,3E20.12)') step, t,  norm_U, norm_V
             !  end if
        end if
    
        if (flag=="f" .or. flag=="F") then ! finish ====================
        end if   

        end select
        abort_flag == .false.
    end subroutine solution_out  
    
        
end module tssm_gray_scott

program test_gray_scott
    use tssm_gray_scott
    implicit none
    
    real(prec) :: t0, tend, dt

    type(gray_scott) :: method 
    type(wf_gray_scott) :: psi    

    call initialize_tssm
 
    
    method = gray_scott(nx = 256,         &
                    xmin = 0.0_prec,  &
                    xmax = 2.5_prec,  &
                    ny = 256,         &
                    ymin = 0.0_prec,  &
                    ymax = 2.5_prec,  &
                    boundary_conditions = periodic) ! Neumann boundary conditions

    psi =  wf_gray_scott(method)
    
    call psi%rset(initial_value)

    dt = 0.1_prec
    t0 = 0.0_prec
    tend = 0.1_prec


    call psi%run(dt = dt, &
                 t0 = t0, &
                 steps=nint((tend-t0)/dt), &
                 splitting_scheme = (/ 0.5_prec, 1.0_prec, 0.5_prec /), &
                 start_with_B = .false., &
                 solution_out = solution_out)

    call finalize_tssm
end program test_gray_scott
