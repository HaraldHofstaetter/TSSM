#define S0(x,y)  x ## _ ## y ## D 
#define S1(x,y) S0(x,y)
#define S(x) S1(x,_DIM_)

#define _WAVE_FUNCTION_ wave_function
#define _COMPLEX_OR_REAL_ complex

#ifdef _QUADPRECISION_
module S(tssmq_fourier_rotating) 
    use tssmq_base
    use tssmq_fourier_common ! fftw_free
    use S(tssmq_fourier)
#else
module S(tssm_fourier_rotating) 
    use tssm_base
    use tssm_fourier_common ! fftw_free
    use S(tssm_fourier)
#endif
    implicit none

    private
    public :: S(fourier_rotating), S(wf_fourier_rotating)

    type, extends(S(fourier)) :: S(fourier_rotating)
        real(kind=prec) :: Omega
        real(kind=prec), pointer :: eigenvalues_d1(:) => null()
        real(kind=prec), pointer :: eigenvalues_d2(:) => null()
    contains    
        procedure :: finalize => finalize_method
    end type S(fourier_rotating)

    interface  S(fourier_rotating) ! constructor
        module procedure new_method
    end interface S(fourier_rotating)


    type, extends(S(wf_fourier)) :: S(wf_fourier_rotating)
    contains    
        procedure :: clone
        procedure :: propagate_A
        procedure :: propagate_A_derivative
        procedure :: add_apply_A
        procedure :: propagate_C
        procedure :: propagate_C_derivative
        procedure :: add_apply_C
    
!        procedure :: finalize => finalize_wf
    end type S(wf_fourier_rotating)

    interface S(wf_fourier_rotating) ! constructor
        module procedure new_wf
    end interface S(wf_fourier_rotating)

contains 

#if(_DIM_==2)
    function new_method(nx, xmin, xmax, ny, ymin, ymax, Omega) result(this)
#elif(_DIM_==3)
    function new_method(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax, Omega) result(this)
#endif
        type(S(fourier_rotating)) :: this
        real(kind=prec),  intent(in) :: xmin 
        real(kind=prec),  intent(in) :: xmax
        integer, intent(in) :: nx
#if(_DIM_>=2)
        real(kind=prec),  intent(in) :: ymin 
        real(kind=prec),  intent(in) :: ymax
        integer, intent(in) :: ny
#endif
#if(_DIM_>=3)
        real(kind=prec),  intent(in) :: zmin 
        real(kind=prec),  intent(in) :: zmax
        integer, intent(in) :: nz
#endif
        real(kind=prec), intent(in) :: Omega

        this%Omega = Omega
#if(_DIM_==2)
        this%S(fourier) = S(fourier)(nx, xmin, xmax, ny, ymin, ymax, boundary_conditions=PERIODIC)
#elif(_DIM_==3)
        this%S(fourier) = S(fourier)(nx, xmin, xmax, ny, ymin, ymax, nz, zmin, zmax, boundary_conditions=PERIODIC)
#endif
        allocate( this%eigenvalues_d1(this%nf1min:this%nf1max ) )
        allocate( this%eigenvalues_d2(this%nf2min:this%nf2max ) )

        this%eigenvalues_d1 = sqrt(-this%eigenvalues1)
        this%eigenvalues_d2 = sqrt(-this%eigenvalues2)
    end function new_method


    subroutine finalize_method(this)
        class(S(fourier_rotating)), intent(inout) :: this

        deallocate( this%eigenvalues_d1 )
        deallocate( this%eigenvalues_d2 )
        call this%S(fourier)%finalize()
    end subroutine finalize_method


    function new_wf(m, coefficient) result(this)
        type(S(wf_fourier_rotating)) :: this
        class(S(fourier_rotating)), target, intent(inout) :: m
        _COMPLEX_OR_REAL_(kind=prec), optional, intent(in) :: coefficient

        this%S(wf_fourier) = S(wf_fourier)(m, coefficient=coefficient)
    end function new_wf


    function clone(this) 
        class(S(wf_fourier_rotating)), intent(inout) :: this
        class(_WAVE_FUNCTION_), pointer :: clone
        type(S(wf_fourier_rotating)), pointer :: p

        select type (m=>this%m); class is (S(fourier_rotating))
        allocate( p )
        p = S(wf_fourier_rotating)(m)
        clone => p
        end select
    end function clone


    subroutine propagate_A(this, dt) ! derivatives w.r.t. x
        implicit none
        class(S(wf_fourier_rotating)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt

#ifdef _OPENMP
#if(_DIM_==2)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf(:,:)
#elif(_DIM_==3)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf(:,:,:)
#endif
        real(kind=prec), pointer :: ev(:)
        real(kind=prec), pointer :: nodes(:)
        integer :: j
#endif        
        if (dt==0.0_prec) then
            return
        end if    

        if (this%m%propagate_time_together_with_A) then
            call this%propagate_time(0.5_prec*dt) 
            ! other half of time propgated by propagate_C ...
        end if

        call this%to_frequency_space_x
        call this%to_real_space_y
#if(_DIM_==3)            
        call this%to_frequency_space_z
#endif        
        select type (m=>this%m); class is (S(fourier_rotating))
#if(_DIM_==2)            
#ifndef _OPENMP
        this%uf = spread(exp((dt*this%coefficient)*m%eigenvalues1),2, m%nf2max-m%nf2min+1) * this%uf
        this%uf = exp(spread((cmplx(0.0_prec, -1.0_prec, kind=prec)*dt*m%Omega)*m%g%nodes_y, &
                             1, m%nf1max-m%nf1min+1)  &
                     *spread(m%eigenvalues_d1, 2, m%nf2max-m%nf2min+1)) * this%uf
#else 
!$OMP PARALLEL DO PRIVATE(j, uf, nodes) 
        do j=1,n_threads
            uf => this%uf(:,lbound(this%uf,2)+m%jf(j-1):lbound(this%uf,2)+m%jf(j)-1)
            nodes => m%g%nodes_y(lbound(m%g%nodes_y,1)+m%jf(j-1):&
                                      lbound(m%g%nodes_y,1)+m%jf(j)-1)
            uf = spread(exp((dt*this%coefficient)*m%eigenvalues1),2, m%jf(j)-m%jf(j-1)) * uf
            uf = exp(spread((cmplx(0.0_prec, -1.0_prec, kind=prec)*dt*m%Omega)*nodes,1, m%nf1max-m%nf1min+1) & 
                     *spread(m%eigenvalues_d1, 2, m%jf(j)-m%jf(j-1))) * uf
        end do
!$OMP END PARALLEL DO 
#endif
#elif(_DIM_==3)            
#ifndef _OPENMP
        this%uf = spread(spread(exp((dt*this%coefficient)*m%eigenvalues1),2, m%nf2max-m%nf2min+1), &
                         3, m%nf3max-m%nf3min+1) * this%uf
        this%uf = spread(spread(exp((0.5_prec*dt*this%coefficient)*m%eigenvalues3), 1, m%nf1max-m%nf1min+1), &
                                2, m%nf2max-m%nf2min+1) * this%uf ! other half in propagate_C
        this%uf = exp(spread(spread((cmplx(0.0_prec, -1.0_prec, kind=prec)*dt*m%Omega)*m%g%nodes_y, &
                             1, m%nf1max-m%nf1min+1)  &
                     *spread(m%eigenvalues_d1, 2, m%nf2max-m%nf2min+1), &
                         3, m%nf3max-m%nf3min+1)) * this%uf
#else 
!$OMP PARALLEL DO PRIVATE(j, uf, ev) 
        do j=1,n_threads
            uf => this%uf(:,:,lbound(this%uf,3)+m%jf(j-1):lbound(this%uf,3)+m%jf(j)-1)
            ev => m%eigenvalues3(lbound(m%eigenvalues3,1)+m%jf(j-1):&
                                       lbound(m%eigenvalues3,1)+m%jf(j)-1)
            uf = spread(spread(exp((dt*this%coefficient)*m%eigenvalues1), &
                            2, m%nf2max-m%nf2min+1), &
                            3,  m%jf(j)-m%jf(j-1)) * uf
            uf = spread(spread(exp((0.5_prec*dt*this%coefficient)*ev), &
                            1, m%nf1max-m%nf1min+1), &
                            2, m%nf2max-m%nf2min+1) * uf
            uf = exp(spread(spread((cmplx(0.0_prec, -1.0_prec, kind=prec)*dt*m%Omega)*m%g%nodes_y, &
                         1, m%nf1max-m%nf1min+1)  &
                 *spread(m%eigenvalues_d1, 2, m%nf2max-m%nf2min+1), &
                     3,  m%jf(j)-m%jf(j-1))) * uf 
        end do
!$OMP END PARALLEL DO 
#endif
#endif
        end select
    end subroutine propagate_A



    subroutine propagate_C(this, dt) ! derivatives w.r.t. y
        implicit none
        class(S(wf_fourier_rotating)), intent(inout) :: this
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt

#ifdef _OPENMP
#if(_DIM_==2)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf(:,:)
#elif(_DIM_==3)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf(:,:,:)
#endif
        real(kind=prec), pointer :: ev(:)
        real(kind=prec), pointer :: evd(:)
        integer :: j
#endif        
        if (dt==0.0_prec) then
            return
        end if    

        if (this%m%propagate_time_together_with_A) then
            call this%propagate_time(0.5_prec*dt) 
            ! other half of time propgated by propagate_A ...
        end if

        call this%to_real_space_x
        call this%to_frequency_space_y
#if(_DIM_==3)            
        call this%to_frequency_space_z
#endif        
        select type (m=>this%m); class is (S(fourier_rotating))
#if(_DIM_==2)            
#ifndef _OPENMP
        this%uf = spread(exp((dt*this%coefficient)*m%eigenvalues2),1, m%nf1max-m%nf1min+1) * this%uf
        this%uf = exp(spread((cmplx(0.0_prec, +1.0_prec, kind=prec)*dt*m%Omega)*m%g%nodes_x, &
                             2, m%nf2max-m%nf2min+1)  &
                     *spread(m%eigenvalues_d2, 1, m%nf1max-m%nf1min+1)) * this%uf
#else 
!$OMP PARALLEL DO PRIVATE(j, uf, ev, evd) 
        do j=1,n_threads
            uf => this%uf(:,lbound(this%uf,2)+m%jf(j-1):lbound(this%uf,2)+m%jf(j)-1)
            ev => m%eigenvalues2(lbound(m%eigenvalues2,1)+m%jf(j-1):&
                                      lbound(m%eigenvalues2,1)+m%jf(j)-1)
            evd => m%eigenvalues2(lbound(m%eigenvalues_d2,1)+m%jf(j-1):&
                                      lbound(m%eigenvalues_d2,1)+m%jf(j)-1)
            uf = spread(exp((dt*this%coefficient)*ev),1, m%nf1max-m%nf1min+1) * uf
            uf = exp(spread((cmplx(0.0_prec, +1.0_prec, kind=prec)*dt*m%Omega)*m%g%nodes_x,&
                             2,m%jf(j)-m%jf(j-1)) & 
                     *spread(evd, 1, m%nf1max-m%nf1min+1)) * uf
        end do
!$OMP END PARALLEL DO 
#endif
#elif(_DIM_==3)            
#ifndef _OPENMP
        this%uf = spread(spread(exp((dt*this%coefficient)*m%eigenvalues2),1, m%nf1max-m%nf1min+1), &
                         3, m%nf3max-m%nf3min+1) * this%uf
        this%uf = spread(spread(exp((0.5_prec*dt*this%coefficient)*m%eigenvalues3), 1, m%nf1max-m%nf1min+1), &
                                2, m%nf2max-m%nf2min+1) * this%uf ! other half in propagate_A
        this%uf = exp(spread(spread((cmplx(0.0_prec, +1.0_prec, kind=prec)*dt*m%Omega)*m%g%nodes_x, &
                             2, m%nf2max-m%nf2min+1)  &
                     *spread(m%eigenvalues_d2, 1, m%nf1max-m%nf1min+1), &
                         3, m%nf3max-m%nf3min+1)) * this%uf
#else 
!$OMP PARALLEL DO PRIVATE(j, uf, ev) 
        do j=1,n_threads
            uf => this%uf(:,:,lbound(this%uf,3)+m%jf(j-1):lbound(this%uf,3)+m%jf(j)-1)
            ev => m%eigenvalues3(lbound(m%eigenvalues3,1)+m%jf(j-1):&
                                       lbound(m%eigenvalues3,1)+m%jf(j)-1)
            uf = spread(spread(exp((dt*this%coefficient)*m%eigenvalues2), &
                            1, m%nf1max-m%nf1min+1), &
                            3,  m%jf(j)-m%jf(j-1)) * uf
            uf = spread(spread(exp((0.5*dt*this%coefficient)*ev), &
                            1, m%nf1max-m%nf1min+1), &
                            2, m%nf2max-m%nf2min+1) * uf
            uf = exp(spread(spread((cmplx(0.0_prec, +1.0_prec, kind=prec)*dt*m%Omega)*m%g%nodes_x, &
                         2, m%nf2max-m%nf2min+1)  &
                 *spread(m%eigenvalues_d2, 1, m%nf1max-m%nf1min+1), &
                     3,  m%jf(j)-m%jf(j-1))) * uf 
        end do
!$OMP END PARALLEL DO 
#endif
#endif
        end select
    end subroutine propagate_C


    subroutine propagate_A_derivative(this, wf, dt)
        class(S(wf_fourier_rotating)), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: wf 
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt
        logical :: f

        select type (wf)
        class is (S(wf_fourier_rotating))
        if (.not.associated(wf%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    

        f = this%m%propagate_time_together_with_A
        this%m%propagate_time_together_with_A = .false.
        call wf%propagate_A(dt)
        this%m%propagate_time_together_with_A = f

        call this%propagate_A(dt)
        
        class default
           stop "E: wrong spectral method for schroedinger wave function"
        end select
    end subroutine propagate_A_derivative


    subroutine propagate_C_derivative(this, wf, dt)
        class(S(wf_fourier_rotating)), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: wf 
        _COMPLEX_OR_REAL_(kind=prec), intent(in) :: dt
        logical :: f

        select type (wf)
        class is (S(wf_fourier_rotating))
        if (.not.associated(wf%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    

        f = this%m%propagate_time_together_with_A
        this%m%propagate_time_together_with_A = .false.
        call wf%propagate_C(dt)
        this%m%propagate_time_together_with_A = f

        call this%propagate_C(dt)
        
        class default
           stop "E: wrong spectral method for schroedinger wave function"
        end select
    end subroutine propagate_C_derivative
    

    subroutine add_apply_A(this, wf, coefficient)
        class(S(wf_fourier_rotating)), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: wf 
        _COMPLEX_OR_REAL_(kind=prec), intent(in), optional :: coefficient
        _COMPLEX_OR_REAL_(kind=prec) :: C 
        _COMPLEX_OR_REAL_(kind=prec) :: Omega 

#ifdef _OPENMP
#if(_DIM_==2)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf1(:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf2(:,:)
#elif(_DIM_==3)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf1(:,:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf2(:,:,:)
#endif
        real(kind=prec), pointer :: ev(:)
        real(kind=prec), pointer :: nodes(:)
        integer :: j
#endif        

        select type (wf)
        class is (S(wf_fourier_rotating))
        
        if (.not.associated(wf%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    

        select type (m=>this%m); class is (S(fourier_rotating))

        C = this%coefficient
        Omega = m%Omega
        if (present(coefficient)) then
            C = C*coefficient
            Omega = Omega*coefficient
        end if    
        if (C==0) then
            return
        end if    

        if (m%propagate_time_together_with_A) then
            if (present(coefficient)) then
                call wf%propagate_time(0.5_prec*coefficient)
            else
#ifdef _REAL_            
                call wf%propagate_time(0.5_prec)
#else                
                call wf%propagate_time((0.5_prec,0.0_prec))
#endif                
            end if
        end if
        
        call this%to_frequency_space_x
        call this%to_real_space_y
        call wf%to_frequency_space_x
        call wf%to_real_space_y
#if(_DIM_==3)            
        call this%to_frequency_space_z
        call wf%to_frequency_space_z
#endif        

#if(_DIM_==2)            
#ifndef _OPENMP            
        wf%uf = wf%uf + spread(C*m%eigenvalues1, 2, m%nf2max-m%nf2min+1) * this%uf
        wf%uf = wf%uf + spread((cmplx(0.0_prec, -1.0_prec, kind=prec)*Omega)*m%g%nodes_y, &
                             1, m%nf1max-m%nf1min+1)  &
                       *spread(m%eigenvalues_d1, 2, m%nf2max-m%nf2min+1) * this%uf
#else
!$OMP PARALLEL DO PRIVATE(j, uf1, uf2, nodes) 
        do j=1,n_threads
            uf1 => wf%uf(:,lbound(wf%uf,2)+m%jf(j-1):lbound(wf%uf,2)+m%jf(j)-1)
            uf2 => this%uf(:,lbound(this%uf,2)+m%jf(j-1):lbound(this%uf,2)+m%jf(j)-1)
            nodes => m%g%nodes_y(lbound(m%g%nodes_y,1)+m%jf(j-1):&
                                      lbound(m%g%nodes_y,1)+m%jf(j)-1)
            uf1 = uf1 + spread(C*m%eigenvalues1, 2,  m%jf(j)-m%jf(j-1)) * uf2
            uf1 = uf1 + spread((cmplx(0.0_prec, -1.0_prec, kind=prec)*Omega)*nodes,1, m%nf1max-m%nf1min+1) & 
                       *spread(m%eigenvalues_d1, 2, m%jf(j)-m%jf(j-1)) * uf2
            end do
!$OMP END PARALLEL DO 
#endif
#elif(_DIM_==3)            
#ifndef _OPENMP  
        wf%uf = wf%uf + spread(spread(C*m%eigenvalues1, 2, m%nf2max-m%nf2min+1), &
                                3, m%nf3max-m%nf3min+1) * this%uf
        wf%uf = wf%uf + spread(spread(0.5_prec*C*m%eigenvalues3, 1, m%nf1max-m%nf1min+1), &
                                2, m%nf2max-m%nf2min+1) * this%uf
        wf%uf = wf%uf + spread(spread((cmplx(0.0_prec, -1.0_prec, kind=prec)*Omega)*m%g%nodes_y, &
                             1, m%nf1max-m%nf1min+1)  &
                     *spread(m%eigenvalues_d1, 2, m%nf2max-m%nf2min+1), &
                         3, m%nf3max-m%nf3min+1) * this%uf
#else
!$OMP PARALLEL DO PRIVATE(j, uf1, uf2, ev) 
        do j=1,n_threads
            uf1 => wf%uf(:,:,lbound(wf%uf,3)+m%jf(j-1):lbound(wf%uf,3)+m%jf(j)-1)
            uf2 => this%uf(:,:,lbound(this%uf,3)+m%jf(j-1):lbound(this%uf,3)+m%jf(j)-1)
            ev => m%eigenvalues3(lbound(m%eigenvalues3,1)+m%jf(j-1):&
                                       lbound(m%eigenvalues3,1)+m%jf(j)-1)
            uf1 = uf1 + C*spread(spread(C*m%eigenvalues1, &
                            2, m%nf2max-m%nf2min+1), &
                            3, m%jf(j)-m%jf(j-1)) * uf2
            uf1 = uf1 + C*spread(spread(C*ev, &
                            1, m%nf1max-m%nf1min+1), &
                            2, m%nf2max-m%nf2min+1) * uf2
            uf1 = uf1 + spread(spread((cmplx(0.0_prec, -1.0_prec, kind=prec)*Omega)*m%g%nodes_y, &
                         1, m%nf1max-m%nf1min+1)  &
                 *spread(m%eigenvalues_d1, 2, m%nf2max-m%nf2min+1), &
                     3,  m%jf(j)-m%jf(j-1)) * uf2
            end do
!$OMP END PARALLEL DO 
#endif

#endif
        end select
        end select
    end subroutine add_apply_A


    subroutine add_apply_C(this, wf, coefficient)
        class(S(wf_fourier_rotating)), intent(inout) :: this
        class(_WAVE_FUNCTION_), intent(inout) :: wf 
        _COMPLEX_OR_REAL_(kind=prec), intent(in), optional :: coefficient
        _COMPLEX_OR_REAL_(kind=prec) :: C 
        _COMPLEX_OR_REAL_(kind=prec) :: Omega 

#ifdef _OPENMP
#if(_DIM_==2)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf1(:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf2(:,:)
#elif(_DIM_==3)            
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf1(:,:,:)
        _COMPLEX_OR_REAL_(kind=prec), pointer :: uf2(:,:,:)
#endif
        real(kind=prec), pointer :: ev(:)
        real(kind=prec), pointer :: evd(:)
        integer :: j
#endif        

        select type (wf)
        class is (S(wf_fourier_rotating))
        
        if (.not.associated(wf%m,this%m)) then
            stop "E: wave functions not belonging to the same method"
        end if    

        select type (m=>this%m); class is (S(fourier_rotating))

        C = this%coefficient
        Omega = m%Omega
        if (present(coefficient)) then
            C = C*coefficient
            Omega = Omega*coefficient
        end if    
        if (C==0) then
            return
        end if    

        if (m%propagate_time_together_with_A) then
            if (present(coefficient)) then
                call wf%propagate_time(0.5_prec*coefficient)
            else
#ifdef _REAL_            
                call wf%propagate_time(0.5_prec)
#else                
                call wf%propagate_time((0.5_prec,0.0_prec))
#endif                
            end if
        end if

        call this%to_real_space_x
        call this%to_frequency_space_y
        call wf%to_real_space_x
        call wf%to_frequency_space_y
#if(_DIM_==3)            
        call this%to_frequency_space_z
        call wf%to_frequency_space_z
#endif                

#if(_DIM_==2)            
#ifndef _OPENMP            
        wf%uf = wf%uf + spread(C*m%eigenvalues2,1, m%nf1max-m%nf1min+1) * this%uf
        wf%uf = wf%uf + spread((cmplx(0.0_prec, +1.0_prec, kind=prec)*Omega)*m%g%nodes_x, &
                             2, m%nf2max-m%nf2min+1)  &
                      *spread(m%eigenvalues_d2, 1, m%nf1max-m%nf1min+1) * this%uf
#else
!$OMP PARALLEL DO PRIVATE(j, uf1, uf2, ev, evd) 
        do j=1,n_threads
            uf1 => wf%uf(:,lbound(wf%uf,2)+m%jf(j-1):lbound(wf%uf,2)+m%jf(j)-1)
            uf2 => this%uf(:,lbound(this%uf,2)+m%jf(j-1):lbound(this%uf,2)+m%jf(j)-1)
            ev => m%eigenvalues2(lbound(m%eigenvalues2,1)+m%jf(j-1):&
                                      lbound(m%eigenvalues2,1)+m%jf(j)-1)
            evd => m%eigenvalues2(lbound(m%eigenvalues_d2,1)+m%jf(j-1):&
                                      lbound(m%eigenvalues_d2,1)+m%jf(j)-1)
            uf1 = uf1 + spread(C*ev, 1, m%nf1max-m%nf1min+1) * uf2
            uf1 = uf1 + spread((cmplx(0.0_prec, +1.0_prec, kind=prec)*Omega)*m%g%nodes_x, &
                             2, m%jf(j)-m%jf(j-1))  &
                      *spread(evd, 1, m%nf1max-m%nf1min+1) * uf2
            end do
!$OMP END PARALLEL DO 
#endif
#elif(_DIM_==3)  
#ifndef _OPENMP  
        wf%uf = wf%uf + spread(spread(C*m%eigenvalues2, 1, m%nf1max-m%nf1min+1), &
                                3, m%nf3max-m%nf3min+1) * this%uf
        wf%uf = wf%uf + spread(spread(0.5_prec*m%eigenvalues3, 1, m%nf1max-m%nf1min+1), &
                                2, m%nf2max-m%nf2min+1) * this%uf
        wf%uf = wf%uf + spread(spread((cmplx(0.0_prec, +1.0_prec, kind=prec)*Omega)*m%g%nodes_x, &
                             2, m%nf2max-m%nf2min+1)  &
                     *spread(m%eigenvalues_d2, 1, m%nf1max-m%nf1min+1), &
                         3, m%nf3max-m%nf3min+1) * this%uf
#else
!$OMP PARALLEL DO PRIVATE(j, uf1, uf2, ev) 
        do j=1,n_threads
            uf1 => wf%uf(:,:,lbound(wf%uf,3)+m%jf(j-1):lbound(wf%uf,3)+m%jf(j)-1)
            uf2 => this%uf(:,:,lbound(this%uf,3)+m%jf(j-1):lbound(this%uf,3)+m%jf(j)-1)
            ev => m%eigenvalues3(lbound(m%eigenvalues3,1)+m%jf(j-1):&
                                       lbound(m%eigenvalues3,1)+m%jf(j)-1)
            uf1 = uf1 + spread(spread(C*m%eigenvalues2, &
                            1, m%nf1max-m%nf1min+1), &
                            3,  m%jf(j)-m%jf(j-1)) * uf2
            uf1 = uf1 + spread(spread((0.5*C)*ev, &
                            1, m%nf1max-m%nf1min+1), &
                            2, m%nf2max-m%nf2min+1) * uf2
            uf1 = uf1 + spread(spread((cmplx(0.0_prec, +1.0_prec, kind=prec)*Omega)*m%g%nodes_x, &
                         2, m%nf2max-m%nf2min+1)  &
                 *spread(m%eigenvalues_d2, 1, m%nf1max-m%nf1min+1), &
                     3,  m%jf(j)-m%jf(j-1)) * uf2 
            end do
!$OMP END PARALLEL DO 
#endif
#endif
        end select
        end select
    end subroutine add_apply_C

    

#ifdef _QUADPRECISION_
end module S(tssmq_fourier_rotating) 
#else
end module S(tssm_fourier_rotating) 
#endif
