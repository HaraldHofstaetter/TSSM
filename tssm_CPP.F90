#ifdef _REAL_
  #define _WAVE_FUNCTION_ real_wave_function
  #define _COMPLEX_OR_REAL_ real
  #define S(x)  x ## _real
#else
  #define _WAVE_FUNCTION_ wave_function
  #define _COMPLEX_OR_REAL_ complex 
  #define S(x)  x ## _complex
#endif  

module S(tssm)
     use tssm_common
     implicit none

     type, abstract :: _WAVE_FUNCTION_
        logical :: is_real_space = .true.
     contains 
        procedure :: propagate_A => S(propagate_A_wf)
        procedure :: propagate_B => S(propagate_B_wf)
        procedure :: propagate_C => S(propagate_C_wf)
        procedure :: step => S(step_wf)
        procedure :: run => S(run_wf)
        procedure :: run_adaptive => S(run_adaptive_wf)
        procedure :: print_orders => S(print_orders_wf)
        procedure :: print_local_orders => S(print_local_orders_wf)
!        procedure(S(proagate_AB_interface)), deferred :: propagate_A 
!        procedure(S(proagate_AB_interface)), deferred :: propagate_B
!        procedure(S(to_rf_space_interface)), deferred :: to_real_space 
!        procedure(S(to_rf_space_interface)), deferred :: to_frequency_space 
        procedure(S(save_load_interface)), deferred :: save
        procedure(S(save_load_interface)), deferred :: load
        procedure(S(norm2_interface)), deferred :: norm2
        procedure(S(distance_interface)), deferred :: distance
        procedure(S(normalize_interface)), deferred :: normalize
        procedure(S(axpy_interface)), deferred :: axpy
        procedure(S(scale_interface)), deferred :: scale
!        procedure(S(set_interface)), deferred :: set
        procedure(S(clone_interface)), deferred :: clone 
        procedure(S(copy_interface)), deferred :: copy
        procedure(S(finalize_interface)), deferred :: finalize 
     end type _WAVE_FUNCTION_

     abstract interface
        subroutine S(to_rf_space_interface)(this)
            import _WAVE_FUNCTION_ 
            class(_WAVE_FUNCTION_), intent(inout) :: this
        end subroutine S(to_rf_space_interface)
        subroutine S(proagate_AB_interface)(this, dt)
            import _WAVE_FUNCTION_, prec 
            class(_WAVE_FUNCTION_), intent(inout) :: this
            _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt
        end subroutine S(proagate_AB_interface)
        subroutine S(save_load_interface)(this, filename)
            import  _WAVE_FUNCTION_ 
            class(_WAVE_FUNCTION_), intent(inout) :: this
            character(len=*), intent(in) :: filename
        end subroutine S(save_load_interface)
        subroutine S(set_interface)(this, f)
            import _WAVE_FUNCTION_, prec 
            class(_WAVE_FUNCTION_), intent(inout) :: this
            _COMPLEX_OR_REAL_(kind=prec), external :: f
        end subroutine S(set_interface)
        function S(norm2_interface)(this)
            import _WAVE_FUNCTION_, prec
            class(_WAVE_FUNCTION_), intent(inout) :: this
            real(kind=prec) :: S(norm2_interface)
        end function S(norm2_interface)
        subroutine S(normalize_interface)(this, norm)
            import _WAVE_FUNCTION_, prec
            class(_WAVE_FUNCTION_), intent(inout) :: this
            real(kind=prec), intent(out) , optional :: norm
        end subroutine S(normalize_interface)
        subroutine S(scale_interface)(this, factor)
            import _WAVE_FUNCTION_, prec
            class(_WAVE_FUNCTION_), intent(inout) :: this
            _COMPLEX_OR_REAL_(kind=prec), intent(in) :: factor
        end subroutine S(scale_interface)
        subroutine S(axpy_interface)(this, other, factor)
            import _WAVE_FUNCTION_, prec
            class(_WAVE_FUNCTION_), intent(inout) :: this
            class(_WAVE_FUNCTION_), intent(inout) :: other
            _COMPLEX_OR_REAL_(kind=prec), intent(in) :: factor
        end subroutine S(axpy_interface)
        function S(distance_interface)(this, wf)
            import  _WAVE_FUNCTION_, prec 
            class(_WAVE_FUNCTION_), intent(inout) :: this
            class(_WAVE_FUNCTION_), intent(inout) :: wf
            real(kind=prec) :: S(distance_interface)
        end function S(distance_interface)
        function S(clone_interface)(this) result(clone)
            import  _WAVE_FUNCTION_ 
            class(_WAVE_FUNCTION_), intent(inout) :: this
            class(_WAVE_FUNCTION_), pointer :: clone
        end function S(clone_interface)
        subroutine S(copy_interface)(this, source)
            import  _WAVE_FUNCTION_ 
            class(_WAVE_FUNCTION_), intent(inout) :: this
            class(_WAVE_FUNCTION_), intent(inout) :: source 
        end subroutine S(copy_interface)
        subroutine S(finalize_interface)(this)
            import  _WAVE_FUNCTION_ 
            class(_WAVE_FUNCTION_), intent(inout) :: this
        end subroutine S(finalize_interface)
    end interface

    interface
        subroutine S(solution_out_interface)(psi, t, dt, flag, abort_flag)
            import :: _WAVE_FUNCTION_, prec
            class(_WAVE_FUNCTION_), intent(inout) :: psi 
            real(kind=prec), intent(in) :: t
            real(kind=prec), intent(in) :: dt
            character, intent(in) :: flag
            logical, intent(out) :: abort_flag
        end subroutine S(solution_out_interface)

        subroutine S(get_reference_solution_interface)(psi, t)
            import :: _WAVE_FUNCTION_, prec
            class(_WAVE_FUNCTION_), intent(inout) :: psi 
            real(kind=prec), intent(in) :: t
        end subroutine S(get_reference_solution_interface)
    end interface

contains 

    subroutine S(propagate_A_wf)(this, dt)
        class(_WAVE_FUNCTION_), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt
        ! dummy routine, do nothing
    end subroutine S(propagate_A_wf)

    subroutine S(propagate_B_wf)(this, dt)
        class(_WAVE_FUNCTION_), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt
        ! dummy routine, do nothing
    end subroutine S(propagate_B_wf)

   subroutine S(propagate_C_wf)(this, dt)
        class(_WAVE_FUNCTION_), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt
        ! dummy routine, do nothing
    end subroutine S(propagate_C_wf)


#ifdef _REAL_
    subroutine S(step_wf)(this, dt, t0, steps, splitting_scheme, operator_sequence, solution_out) 
#else    
    subroutine S(step_wf)(this, dt, t0, steps, splitting_scheme, complex_splitting_scheme, operator_sequence, solution_out) 
#endif    
        class(_WAVE_FUNCTION_), intent(inout) :: this
        real(kind=prec), intent(in) :: dt
        real(kind=prec), intent(in) :: t0 
        integer, intent(in), optional :: steps
        character(len=*), intent(in), optional :: operator_sequence
        real(kind=prec), intent(in), target, optional :: splitting_scheme(:)
#ifndef _REAL_        
        complex(kind=prec), intent(in), target, optional :: complex_splitting_scheme(:)
#endif        
        procedure(S(solution_out_interface)), optional :: solution_out

        integer :: nsteps = 1
        integer :: op, ops
        integer :: split_steps
        integer :: step, split_step
        real(kind=prec) :: t
        logical :: abort_flag
        procedure(S(solution_out_interface)), pointer :: sol_out => null()
        character(len=10) :: seq
#ifdef _REAL_
        real(kind=prec), pointer :: scheme(:)
        real(kind=prec), target :: default_scheme(3) =(/ 0.5_prec, 1.0_prec, 0.5_prec /)
#else
        complex(kind=prec), pointer :: scheme(:)
        complex(kind=prec), target :: default_scheme(3) = & 
                           (/ (0.5_prec,0.0_prec), (1.0_prec,0.0_prec), (0.5_prec,0.0_prec) /)
        integer :: i
#endif 
        if(present(splitting_scheme)) then
#ifdef _REAL_
            scheme => splitting_scheme
#else
            allocate( scheme(ubound(splitting_scheme,1)) )
            !scheme = (/ (cmplx(splitting_scheme(i), 0.0_prec, kind=prec), i=1,ubound(splitting_scheme,1)) /)
            scheme = cmplx(splitting_scheme, 0.0_prec, kind=prec)
        elseif(present(complex_splitting_scheme)) then
            scheme => complex_splitting_scheme
#endif             
        else    
            scheme => default_scheme
        end if    


        ! TODO argument checking

        if (present(steps)) then
            nsteps = steps
        end if

        if (present(operator_sequence)) then
           seq = operator_sequence
           ops = len(trim(operator_sequence))
        else
           seq = "AB"
           ops = 2
        end if     

        if (present(solution_out)) then
           sol_out => solution_out
        end if   

        split_steps = ubound(scheme,1)

        t = t0

        do step=1, nsteps
            op = 1 
            do split_step=1, split_steps
                select case(seq(op:op))
                case("A")
                    call this%propagate_A(dt*scheme(split_step))
                case("B")
                    call this%propagate_B(dt*scheme(split_step))
                case("C")
                    call this%propagate_C(dt*scheme(split_step))
                end select 
                op = op + 1
                if (op>ops) op = 1
            end do
            t = t + dt
            if (associated(sol_out)) then
                call sol_out(this, t, dt, "execute", abort_flag)
                if (abort_flag) exit
            end if    

        end do
#ifndef _REAL_
        if (present(splitting_scheme)) then
            deallocate( scheme ) 
        end if    
#endif
    end subroutine S(step_wf)


    subroutine S(run_wf)(this, dt, times, splitting_scheme, &
#ifndef _REAL_
                 complex_splitting_scheme, &
#endif
                 operator_sequence, solution_out) 
        class(_WAVE_FUNCTION_), intent(inout) :: this
        real(kind=prec), intent(inout) :: dt
        real(kind=prec), intent(in) :: times(:) 
        character(len=*), intent(in), optional :: operator_sequence
        real(kind=prec), intent(in) :: splitting_scheme(:)
#ifndef _REAL_
        complex(kind=prec), intent(in), optional :: complex_splitting_scheme(:)
#endif
        procedure(S(solution_out_interface)), optional :: solution_out

        logical :: abort_flag
        integer :: step, steps
        real(kind=prec) :: t0, tend

        if (present(solution_out)) then
            call solution_out(this, times(1), dt, "init", abort_flag)
        end if

        do step = 2, ubound(times,1)
            t0 = times(step-1)
            tend = times(step)
            steps = ceiling(abs(tend-t0)/dt)
            dt = (tend-t0)/real(steps, kind=prec)
            call S(step_wf)(this, dt, t0, steps, splitting_scheme, &
#ifndef _REAL_
                 complex_splitting_scheme, &
#endif                  
                 operator_sequence, solution_out) 

            if (present(solution_out)) then
                call solution_out(this, tend, dt, "checkpoint", abort_flag)
                if (abort_flag) exit
            end if
        end do

        if (present(solution_out)) then
            call solution_out(this, times(ubound(times,1)), dt, "finish", abort_flag)
        end if    

    end subroutine S(run_wf)



    subroutine S(run_adaptive_wf)(this, dt, times, tol, &
                 splitting_scheme, controller_scheme, order, operator_sequence, solution_out) 
        class(_WAVE_FUNCTION_), intent(inout) :: this
        real(kind=prec), intent(inout) :: dt
        real(kind=prec), intent(in) :: times(:) 
        real(kind=prec), intent(in) :: tol
        character(len=*), intent(in), optional :: operator_sequence
        real(kind=prec), intent(in) :: splitting_scheme(:)
        real(kind=prec), intent(in), optional :: controller_scheme(:)
        integer, intent(in)  :: order
        procedure(S(solution_out_interface)), optional :: solution_out


        logical :: abort_flag
        integer :: step
        real(kind=prec) :: t0, tend

        if (present(solution_out)) then
            call solution_out(this, times(1), dt, "init", abort_flag)
        end if

        do step = 2, ubound(times,1)
            t0 = times(step-1)
            tend = times(step)
            if (tend<t0) then
                dt = -abs(dt)
            else
                dt = abs(dt)
            endif    
            call S(run_adaptive_wf_0)(this, dt, t0, tend, tol, &
                 splitting_scheme, controller_scheme, order, operator_sequence, solution_out)
            if (present(solution_out)) then
                call solution_out(this, tend, dt, "checkpoint", abort_flag)
                if (abort_flag) exit
            end if
        end do

        if (present(solution_out)) then
            call solution_out(this, times(ubound(times,1)), dt, "finish", abort_flag)
        end if    

    end subroutine S(run_adaptive_wf)

   

    subroutine S(run_adaptive_wf_0)(this, dt, t0, tend, tol, &
                 splitting_scheme, controller_scheme, order, operator_sequence, solution_out) 
        class(_WAVE_FUNCTION_), intent(inout) :: this
        real(kind=prec), intent(inout) :: dt
        real(kind=prec), intent(in) :: t0
        real(kind=prec), intent(in) :: tend
        real(kind=prec), intent(in) :: tol
        character(len=*), intent(in), optional :: operator_sequence
        real(kind=prec), intent(in) :: splitting_scheme(:)
        real(kind=prec), intent(in), optional :: controller_scheme(:)
        integer, intent(in)  :: order
        procedure(S(solution_out_interface)), optional :: solution_out

        !real(kind=prec), parameter :: fac = 0.8_prec
        real(kind=prec), parameter :: facmin = 0.25_prec
        real(kind=prec), parameter :: facmax = 4.0_prec
        real(kind=prec), parameter :: fac = 0.9_prec
        !real(kind=prec), parameter :: facmin = 0.5_prec
        !real(kind=prec), parameter :: facmax = 2.0_prec

        integer :: op, ops, op0_AS
        character(len=10) :: seq
        integer :: split_steps, s, split_steps_AS, first_step_A
        integer :: step, split_step, first_step_AS
        real(kind=prec) :: t, dt_save, err, sign_dt, t_out
        logical :: abort_flag
        character :: controller_type

        class(_WAVE_FUNCTION_), pointer :: psi_AS, psi0
        ! TODO argument checking

        if (present(operator_sequence)) then
           seq = operator_sequence
           ops = len(trim(operator_sequence))
        else
           seq = "AB"
           ops = 2
        end if     

        split_steps = ubound(splitting_scheme,1)

        if (present(controller_scheme)) then
            controller_type = 'e' ! embedded pairs
            split_steps_AS = ubound(controller_scheme,1)
            s = min(split_steps, split_steps_AS)
            first_step_AS = maxloc(merge(1,0,splitting_scheme(1:s)/=controller_scheme(1:s)),1) 
            if (first_step_AS==1.and.splitting_scheme(1)==controller_scheme(1)) then
               first_step_AS = s+1
            endif    
        else    
            controller_type = 'p' ! palindromic pairs
            !TODO check if scheme is indeed palindromic and p is odd
            ! If not then fall back to defect based error estimate 
        end if

        psi0 => this%clone()
        psi_AS => this%clone()

        t = t0 

        if (dt<0) then
            sign_dt = -1.0_prec
        else    
            sign_dt = 1.0_prec
        end if    
       
        do while(sign_dt*(tend-t)>0.0_prec)
            call psi0%copy(this)
            do
                dt_save = dt
                if (dt>0) then
                    dt = min(dt, tend-t)
                else    
                    dt = max(dt, tend-t)
                end if
                op = 1  
                if (controller_type=='p') then ! palindromic pairs
                    call psi_AS%copy(this)
                end if    
                do split_step=1, split_steps
                    if (controller_type=='e'.and.split_step==first_step_AS) then
                         op0_AS = op
                         call psi_AS%copy(this)
                    end if
#ifdef _REAL_                    
                    select case(seq(op:op))
                    case("A")
                        call this%propagate_A(dt*splitting_scheme(split_step))
                    case("B")
                        call this%propagate_B(dt*splitting_scheme(split_step))
                    case("C")
                        call this%propagate_C(dt*splitting_scheme(split_step))
                    end select
#else                    
                    select case(seq(op:op))
                    case("A")
                        call this%propagate_A(cmplx(dt*splitting_scheme(split_step),0.0_prec,prec))
                    case("B")
                        call this%propagate_B(cmplx(dt*splitting_scheme(split_step),0.0_prec,prec))
                    case("C")
                        call this%propagate_C(cmplx(dt*splitting_scheme(split_step),0.0_prec,prec))
                    end select    
#endif                    
                    op = op + 1
                    if (op>ops) op = 1
                end do
                select case(controller_type)
                case("e") ! embedded pairs
                  op = op0_AS
                  do split_step=first_step_AS, split_steps_AS
#ifdef _REAL_                    
                    select case(seq(op:op))
                    case("A")
                        call psi_AS%propagate_A(dt*controller_scheme(split_step))
                    case("B")
                        call psi_AS%propagate_B(dt*controller_scheme(split_step))
                    case("C")
                        call psi_AS%propagate_C(dt*controller_scheme(split_step))
                    end select
#else
                    select case(seq(op:op))
                    case("A")
                        call psi_AS%propagate_A(cmplx(dt*controller_scheme(split_step),0.0_prec,prec))
                    case("B")
                        call psi_AS%propagate_B(cmplx(dt*controller_scheme(split_step),0.0_prec,prec))
                    case("C")
                        call psi_AS%propagate_C(cmplx(dt*controller_scheme(split_step),0.0_prec,prec))
                    end select
#endif
                    op = op + 1
                    if (op>ops) op = 1
                  end do    
                  err = this%distance(psi_AS)/tol
                case("p") ! palindromic pairs
                  op = ops  
                  do split_step=1, split_steps
#ifdef _REAL_                    
                    select case(seq(op:op))
                    case("A")
                        call psi_AS%propagate_A(dt*splitting_scheme(split_step))
                    case("B")
                        call psi_AS%propagate_B(dt*splitting_scheme(split_step))
                    case("C")
                        call psi_AS%propagate_C(dt*splitting_scheme(split_step))
                    end select
#else                    
                    select case(seq(op:op))
                    case("A")
                        call psi_AS%propagate_A(cmplx(dt*splitting_scheme(split_step),0.0_prec,prec))
                    case("B")
                        call psi_AS%propagate_B(cmplx(dt*splitting_scheme(split_step),0.0_prec,prec))
                    case("C")
                        call psi_AS%propagate_C(cmplx(dt*splitting_scheme(split_step),0.0_prec,prec))
                    end select    
#endif                    
                    op = op - 1
                    if (op<1) op = ops
                  end do
                  err = 0.5_prec*this%distance(psi_AS)/tol
                end select


                dt = dt*min(facmax, max(facmin, fac*(1.0_prec/err)**(1.0_prec/(real(order+1,prec)))))
                if (err<=1.0_prec) then 
                    !if (this_proc==0) then
                    !    write (*, '(A,E17.8,A,E17.8,A,E17.8,A)') &
                    !      't=', t, '  err=', err, '  dt=', dt, '  accepted...'
                    !end if      
                    exit
                end if


                call this%copy(psi0)
                if (this_proc==0) then
                    write (*, '(A,E17.8,A,E17.8,A,E17.8,A)')  &
                           '# t=', t, '  err=', err, '  dt=', dt, '  rejected...'
                end if           

            end do  

            if (present(solution_out)) then
                if (dt>0) then
                    t_out = min(t + dt, tend)
                else    
                    t_out = max(t + dt, tend)
                end if    

                call solution_out(this, t_out, dt, "execute", abort_flag)
                if (abort_flag) exit
            end if 

            t = t + dt_save
        end do

        dt = dt_save ! on output dt contains last predicted stepsize

        call psi0%finalize
        call psi_AS%finalize

    end subroutine S(run_adaptive_wf_0)



#ifdef _REAL_
    subroutine S(print_orders_wf)(this, reference_solution, t0, tend, rows, dt, splitting_scheme, operator_sequence) 
#else
    subroutine S(print_orders_wf)(this, reference_solution, t0, tend, rows, dt, &
                                  splitting_scheme, complex_splitting_scheme, operator_sequence) 
#endif
        class(_WAVE_FUNCTION_), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: reference_solution
        real(kind=prec), intent(in) ::t0
        real(kind=prec), intent(in) ::tend 
        integer, intent(in), optional :: rows 
        real(kind=prec), intent(in) :: dt
        real(kind=prec), intent(in),optional :: splitting_scheme(:)
#ifndef _REAL_        
        complex(kind=prec), intent(in),optional :: complex_splitting_scheme(:)
#endif        
        character(len=*), intent(in), optional :: operator_sequence

        character(len=*), parameter :: filename = "print_order_saved_solution.h5"

        integer :: steps, row
        real(kind=prec) :: dt1, err, err_old, p

        class(_WAVE_FUNCTION_), pointer :: wf_save_initial_value
    
        wf_save_initial_value => this%clone()
        call wf_save_initial_value%copy(this)

        steps = floor((tend-t0)/dt)
        dt1 = dt
        if (this_proc==0) then
            write (*,*) "            dt         err      p"
            write (*,*) "----------------------------------"
        end if    
        do row=1,rows
#ifdef _REAL_
            call this%step(dt1, t0, steps, splitting_scheme, operator_sequence)
#else
            call this%step(dt1, t0, steps, splitting_scheme, complex_splitting_scheme, operator_sequence)
#endif        
            err = this%distance(reference_solution)
            if (this_proc==0) then
                if (row==1) then
                    write (*, '(I3,2E12.3)') row, dt1, err
                else
                    p = log(err_old/err)/log(2.0_prec);
                    write (*, '(I3,2E12.3, F7.2)') row, dt1, err, p
                end if
            end if    
            err_old = err;
            dt1 = dt1/2.0_prec
            steps = 2*steps
            call this%copy(wf_save_initial_value)
        end do
        call wf_save_initial_value%finalize
    end subroutine S(print_orders_wf)

#ifdef _REAL_
    subroutine S(print_local_orders_wf)(this, rows, dt, splitting_scheme, operator_sequence, &
                                        controller_scheme, &
                                        fraunz, reference_scheme, get_reference_solution) 
#else
    subroutine S(print_local_orders_wf)(this, rows, dt, &
                                        splitting_scheme, complex_splitting_scheme, operator_sequence, &
                                        controller_scheme, complex_controller_scheme, &
                                        fraunz, reference_scheme, complex_reference_scheme, get_reference_solution) 
#endif
        class(_WAVE_FUNCTION_), intent(inout) :: this
        integer, intent(in), optional :: rows 
        real(kind=prec), intent(in) :: dt
        real(kind=prec), intent(in), target, optional :: splitting_scheme(:)
        real(kind=prec), intent(in), target, optional :: controller_scheme(:)
        real(kind=prec), intent(in), target, optional :: reference_scheme(:)
#ifndef _REAL_        
        complex(kind=prec), intent(in), target, optional :: complex_splitting_scheme(:)
        complex(kind=prec), intent(in), target, optional :: complex_controller_scheme(:)
        complex(kind=prec), intent(in), target, optional :: complex_reference_scheme(:)
#endif        
        character(len=*), intent(in), optional :: operator_sequence
        integer, intent(in), optional :: fraunz
        procedure(S(get_reference_solution_interface)), optional :: get_reference_solution

        integer :: row, fraunz1, rows1
        real(kind=prec) :: dt1, err, err_old, p
        real(kind=prec) :: err1, err_old1, p1
        real(kind=prec) :: err2, err_old2, p2
        real(kind=prec), pointer :: reference_scheme1(:)
#ifndef _REAL_        
        complex(kind=prec), pointer :: complex_reference_scheme1(:)
#endif        
        class(_WAVE_FUNCTION_), pointer :: reference_solution, wf_save_initial_value
        class(_WAVE_FUNCTION_), pointer :: controller_solution => null()
        logical :: palindromic_controller
        character(len=10) :: reversed_operator_sequence
        integer :: ll, i

        fraunz1 = 10
        if (present(fraunz)) then
            fraunz1 = fraunz
        end if    

        rows1 = 10
        if (present(rows)) then
            rows1 = rows
        end if    

        reference_scheme1 => splitting_scheme
        if (present(reference_scheme)) then
           reference_scheme1 => reference_scheme
        end if   
#ifndef _REAL_        
        complex_reference_scheme1 => complex_splitting_scheme
        if (present(complex_reference_scheme)) then
           complex_reference_scheme1 => complex_reference_scheme
        end if   
#endif 
    
        wf_save_initial_value => this%clone()
        reference_solution => this%clone()

#ifdef _REAL_        
        if (present(controller_scheme)) then
#else        
        if (present(controller_scheme).or.present(complex_controller_scheme)) then
#endif        
            controller_solution => this%clone()
        endif

        palindromic_controller =.false. 
        if (present(controller_scheme)) then
            palindromic_controller = all(controller_scheme==palindromic)
        end if    

        if (palindromic_controller) then
           if (present(operator_sequence)) then
               reversed_operator_sequence = operator_sequence
               ll = len(trim(  reversed_operator_sequence))
               forall (i=1:ll) reversed_operator_sequence(i:i) = &
                   reversed_operator_sequence(ll-i+1:ll-i+1)
           else
               reversed_operator_sequence = "BA" ! default operartor_sequence = "AB"
           endif
        endif   

        call wf_save_initial_value%copy(this)        
        call reference_solution%copy(wf_save_initial_value)
        if (associated(controller_solution)) then
            call controller_solution%copy(wf_save_initial_value)
        end if

        dt1 = dt
        if (this_proc==0) then
            if (palindromic_controller) then
                write (*,*) "            dt         err      p          err      p          err      p"
                write (*,*) "--------------------------------------------------------------------------"
            else if (associated(controller_solution)) then
                write (*,*) "            dt         err      p          err      p"
                write (*,*) "------------------------------------------------------"
            else
                write (*,*) "            dt         err      p"
                write (*,*) "----------------------------------"
            end if    
        end if    
        do row=1,rows1
            if (present(get_reference_solution)) then
                call get_reference_solution(reference_solution, dt1)
            else
#ifdef _REAL_
                call reference_solution%step(dt1/real(fraunz1, kind=prec), 0.0_prec, fraunz1, &
                                             splitting_scheme, operator_sequence)
#else           
                if (present(complex_splitting_scheme)) then
                    call reference_solution%step(dt1/real(fraunz1, kind=prec), 0.0_prec, fraunz1, & 
                                             complex_splitting_scheme = complex_reference_scheme1, &
                                             operator_sequence = operator_sequence)
                else 
                    call reference_solution%step(dt1/real(fraunz1, kind=prec), 0.0_prec, fraunz1, & 
                                             splitting_scheme = reference_scheme1, &
                                             operator_sequence = operator_sequence)
                end if
#endif        
            end if

#ifdef _REAL_
            call this%step(dt1, 0.0_prec, 1, splitting_scheme, operator_sequence)
#else
            call this%step(dt1, 0.0_prec, 1, splitting_scheme, complex_splitting_scheme, operator_sequence)
#endif        
            err = this%distance(reference_solution)
            if (palindromic_controller) then
#ifdef _REAL_
                call controller_solution%step(dt1, 0.0_prec, 1, splitting_scheme, reversed_operator_sequence)
                err1 = controller_solution%distance(reference_solution)
                call controller_solution%axpy(this, 1.0_prec)
                call controller_solution%scale(0.5_prec)
#else
                call controller_solution%step(dt1, 0.0_prec, 1, splitting_scheme, complex_splitting_scheme, &
                     reversed_operator_sequence)
                err1 = controller_solution%distance(reference_solution)
                call controller_solution%axpy(this, (1.0_prec, 0.0_prec))
                call controller_solution%scale((0.5_prec, 0.0_prec))
#endif        
                err2 = controller_solution%distance(reference_solution)
            else if (associated(controller_solution)) then
#ifdef _REAL_
                call controller_solution%step(dt1, 0.0_prec, 1, controller_scheme, operator_sequence)
#else
                call controller_solution%step(dt1, 0.0_prec, 1, controller_scheme, complex_controller_scheme, operator_sequence)
#endif        
                err1 = controller_solution%distance(reference_solution)
            end if

            if (this_proc==0) then
                if (row==1) then
                    if (palindromic_controller) then
                        write (*, '(I3,2E12.3,2E20.3)') row, dt1, err, err1, err2
                    else if (associated(controller_solution)) then
                        write (*, '(I3,2E12.3,E20.3)') row, dt1, err, err1
                    else
                        write (*, '(I3,2E12.3)') row, dt1, err
                    end if    
                else
                    p = log(err_old/err)/log(2.0_prec);
                    if (palindromic_controller) then
                        p1 = log(err_old1/err1)/log(2.0_prec);
                        p2 = log(err_old2/err2)/log(2.0_prec);
                        write (*, '(I3,2E12.3, F7.2, E13.3, F7.2, E13.3, F7.2)') row, dt1, err, p, err1, p1, err2, p2
                    else if (associated(controller_solution)) then
                        p1 = log(err_old1/err1)/log(2.0_prec);
                        write (*, '(I3,2E12.3, F7.2, E13.3, F7.2)') row, dt1, err, p, err1, p1
                    else
                        write (*, '(I3,2E12.3, F7.2)') row, dt1, err, p
                    end if
                end if
            end if
            err_old = err;
            call this%copy(wf_save_initial_value)
            call reference_solution%copy(wf_save_initial_value)
            if (associated(controller_solution)) then
                call controller_solution%copy(wf_save_initial_value)
                err_old1 = err1
                if (palindromic_controller) then
                   err_old2 = err2
                end if   
            endif
            dt1 = dt1/2.0_prec
        end do

        call reference_solution%finalize 
        call wf_save_initial_value%finalize

    end subroutine S(print_local_orders_wf)


end module S(tssm)



