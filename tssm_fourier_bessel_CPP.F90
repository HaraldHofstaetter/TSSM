#ifdef _REAL_
#ifdef _QUADPRECISION_
 #define _MODULE_ tssmq_fourier_bessel_real_2d
#else
 #define _MODULE_ tssm_fourier_bessel_real_2d
#endif
 #define _COMPLEX_OR_REAL_ real
 #define _METHOD_ fourier_bessel_real_2d
 #define _WF_ wf_fourier_bessel_real_2d
 #define _BASE_METHOD_ polar_real_2d
 #define _BASE_WF_ wf_polar_real_2d
 #define _METHOD_ROTSYM_ bessel_rotsym_real_1d
 #define _WF_ROTSYM_ wf_bessel_rotsym_real_1d
 #define _BASE_METHOD_ROTSYM_ tensorial_real_1d
 #define _BASE_WF_ROTSYM_ wf_tensorial_real_1d

#else
#ifdef _QUADPRECISION_
 #define _MODULE_ tssmq_fourier_bessel_2d
#else
 #define _MODULE_ tssm_fourier_bessel_2d
#endif
 #define _COMPLEX_OR_REAL_ complex
 #define _METHOD_ fourier_bessel_2d
 #define _WF_ wf_fourier_bessel_2d
 #define _BASE_METHOD_ polar_2d
 #define _BASE_WF_ wf_polar_2d
 #define _METHOD_ROTSYM_ bessel_rotsym_1d
 #define _WF_ROTSYM_ wf_bessel_rotsym_1d
 #define _BASE_METHOD_ROTSYM_ tensorial_1d
 #define _BASE_WF_ROTSYM_ wf_tensorial_1d
#endif


module _MODULE_
#ifdef _QUADPRECISION_
    use tssmq
    use tssmq_grid
    use tssmq_polar
    use tssmq_tensorial
    use tssmq_fourier_common
    use tssmq_fourier_bessel_common
#else
    use tssm
    use tssm_grid
    use tssm_polar
    use tssm_tensorial
    use tssm_fourier_common
    use tssm_fourier_bessel_common
#endif    
    implicit none

    private
    public ::  _METHOD_, _WF_, _METHOD_ROTSYM_, _WF_ROTSYM_

    type, extends(_BASE_METHOD_) :: _METHOD_ 
        real(kind=prec) :: rmax = 1.0_prec
        integer :: boundary_conditions = dirichlet 
        integer :: quadrature_formula = lobatto 
        real(kind=prec), allocatable :: normalization_factors(:,:)
    contains    
        procedure :: finalize => finalize_method 
    end type _METHOD_

    interface _METHOD_ ! constructor
        module procedure new_method
    end interface _METHOD_

    type, extends(_BASE_WF_) ::  _WF_ 
    contains
        procedure :: save
        procedure :: load
        procedure :: evaluate
    end type _WF_ 

    interface _WF_ ! constructor
        module procedure new_wf 
    end interface _WF_

!!!!!!

    type, extends(_BASE_METHOD_ROTSYM_) :: _METHOD_ROTSYM_ 
        real(kind=prec) :: rmax = 1.0_prec
        integer :: boundary_conditions = dirichlet 
!        integer :: quadrature_formula = lobatto 
        real(kind=prec), allocatable :: normalization_factors(:)
    contains    
        procedure :: finalize => finalize_method_rotsym 
    end type _METHOD_ROTSYM_

    interface _METHOD_ROTSYM_ ! constructor
        module procedure new_method_rotsym
    end interface _METHOD_ROTSYM_

    type, extends(_BASE_WF_ROTSYM_) ::  _WF_ROTSYM_ 
    contains
!        procedure :: save => save_rotsym
!        procedure :: load => load_rotsym
!        procedure :: evaluate => evaluate_rotsym
    end type _WF_ROTSYM_ 

    interface _WF_ROTSYM_ ! constructor
        module procedure new_wf_rotsym 
    end interface _WF_ROTSYM_



contains


    function new_method(M, nr, nfr, rmax, boundary_conditions, quadrature_formula) result(this)
        type(_METHOD_) :: this
        integer, intent(in) :: M
        integer, intent(in) :: nr 
        integer, intent(in) :: nfr 
        real(kind=prec), intent(in), optional :: rmax 
        integer, intent(in), optional :: boundary_conditions
        integer, intent(in), optional :: quadrature_formula

        if (present(rmax)) then
            this%rmax = rmax
        end if
        if (present(boundary_conditions)) then
            this%boundary_conditions = boundary_conditions
            if (boundary_conditions==neumann) then
                this%quadrature_formula = radau
            endif     
        end if
        if (present(quadrature_formula)) then
            this%quadrature_formula = quadrature_formula
        endif     

        this%_BASE_METHOD_ = _BASE_METHOD_(M, nr, nfr, .true., .false. ) 
               ! symmetric coefficients
               ! not separated eigenvalues ...

        allocate( this%normalization_factors(this%nfrmin:this%nfrmax, &
                                             this%nfthetamin:this%nfthetamax) ) 

        call fourier_bessel_coeffs(nr, nfr, M, this%g%nodes_r,  this%g%weights_r, this%L, &
                                  this%eigenvalues_r_theta, this%normalization_factors, this%nfthetamin, this%nfthetamax, &
                                  this%boundary_conditions, this%quadrature_formula) 

    end function new_method 

    subroutine finalize_method(this)
        class(_METHOD_), intent(inout) :: this

        call this%_BASE_METHOD_%finalize()
        deallocate( this%normalization_factors )

    end subroutine finalize_method


    function new_wf(m, u, coefficient) result(this)
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_ptr, c_loc
        type(_WF_) :: this
        class(_METHOD_), target, intent(inout) :: m
#ifdef _REAL_       
        real(kind=prec), optional, target, intent(inout) :: u(:,:)
        real(kind=prec), optional, intent(in) :: coefficient
#else            
        complex(kind=prec), optional, target, intent(inout) :: u(:,:)
        complex(kind=prec), optional, intent(in) :: coefficient
#endif
        this%_BASE_WF_ = _BASE_WF_(m, u, coefficient)
    end function new_wf



    subroutine save(this, filename)
#ifdef _NO_HDF5_
        class(_WF_), intent(inout) :: this
        character(len=*), intent(in) :: filename
        print *, "W: save not implemented"
#else
#ifdef _QUADPRECISION_
        use tssmq_hdf5
#else
        use tssm_hdf5
#endif        
        class(_WF_), intent(inout) :: this
        character(len=*), intent(in) :: filename
        !TODO

        call this%_BASE_WF_%save(filename)

        select type (m=>this%m ); class is (_METHOD_)
        call hdf5_write_attributes(filename, (/ "nodes_type" /), &
          (/ "fourier_bessel" /), &
          (/ character(len=0) :: /), (/ integer :: /), &
          (/ "rmax" /), (/ m%rmax /) ) 
        end select
!TODO save nodes !!!
#endif        
    end subroutine save


   subroutine load(this, filename)
#ifdef _NO_HDF5_
        class(_WF_), intent(inout) :: this
        character(len=*), intent(in) :: filename
        print *, "W: load not implemented"
#else
#ifdef _QUADPRECISION_
        use tssmq_hdf5
#else
        use tssm_hdf5
#endif        
        class(_WF_), intent(inout) :: this
        character(len=*), intent(in) :: filename
#ifdef _QUADPRECISION_
        real(kind=prec), parameter :: eps = epsilon(1.0_8)
#else        
        real(kind=prec), parameter :: eps = epsilon(1.0_prec)
#endif
        character(len=20) :: s
        integer :: n, len

       call hdf5_read_string_attribute(filename, "nodes_type", s, len)      
       if (s(1:len)/="fourier_bessel") then
              stop "E: wrong nodes type"
       end if
       select type (m=>this%m ); class is (_METHOD_)
       if (abs(hdf5_read_double_attribute(filename, "rmax")- m%rmax)>eps) then
              stop "E: incompatible grids"
       end if       
       end select

       call this%_BASE_WF_%save(filename)
#endif
    end subroutine load


    function evaluate(this, x, y) result(z)
        class(_WF_), intent(inout) :: this
        real(kind=prec), intent(in) :: x, y
        _COMPLEX_OR_REAL_(kind=prec) :: z

        real(kind=prec) :: r, r2, theta, f, lambda
        complex(kind=prec) :: s
        integer :: j, m, m1, m2

        select type (mm =>this%m ); class is (_METHOD_)

        r2 = x**2 + y**2
        z = 0.0_prec
        if (r2 > mm%rmax) then
            return
        end if

        r = sqrt(r2)
        theta = atan2(y, x)
        call this%to_frequency_space

#ifdef _REAL_                
!$OMP PARALLEL DO PRIVATE(m, m1, m2, j, f, lambda, s) REDUCTION(+:z)
        do m = this%m%nfthetamin, this%m%nfthetamax
            s = 0.0_prec
            do j = this%m%nfrmin, this%m%nfrmax
                lambda = sqrt(this%m%eigenvalues_r_theta(j, m))
                f = mm%normalization_factors(j, m)
                s = s + f*bessel_jn(m, lambda*r) * this%ufc(j, m)
            end do
            if (m/=0 .and. m/=this%m%g%ntheta/2) then
                s = 2*s
            end if
            z = z + real(exp(cmplx(0.0_prec, m*theta, kind=prec)) * s, kind=prec)
        end do
!$OMP END PARALLEL DO        
#else
!$OMP PARALLEL DO PRIVATE(m, m1, m2, j, f, lambda, s) REDUCTION(+:z)
        do m = this%m%nfthetamin, this%m%nfthetamax
            m1 = m
            m2 = m
            if (m1>this%m%g%ntheta/2) then
                m1 = abs(m1-this%m%g%ntheta)
                m2 = -m1
            end if
            s = 0.0_prec
            do j = this%m%nfrmin, this%m%nfrmax
                lambda = sqrt(this%m%eigenvalues_r_theta(j, m1))
                f = mm%normalization_factors(j, m1)
                s = s + f*bessel_jn(m1, lambda*r) * this%uf(j, m)
            end do
            z = z + exp(cmplx(0.0_prec, m2*theta, kind=prec)) * s
        end do
!$OMP END PARALLEL DO        
#endif                
        end select
        
    end function evaluate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    function new_method_rotsym(nr, rmax, boundary_conditions) result(this)
        type(_METHOD_ROTSYM_) :: this
        integer, intent(in) :: nr 
        real(kind=prec), intent(in), optional :: rmax 
        integer, intent(in), optional :: boundary_conditions

        if (present(rmax)) then
            this%rmax = rmax
        end if
        if (present(boundary_conditions)) then
            this%boundary_conditions = boundary_conditions
        end if

        this%_BASE_METHOD_ROTSYM_ = _BASE_METHOD_ROTSYM_(nr) 

        allocate( this%normalization_factors(this%nf1min:this%nf1max) )

        call bessel_rotsym_coeffs(nr, this%g%nodes_x,  this%g%weights_x, this%H1, &
             this%eigenvalues1, this%normalization_factors, this%boundary_conditions) 
    end function new_method_rotsym


     subroutine finalize_method_rotsym(this)
        class(_METHOD_ROTSYM_), intent(inout) :: this

        call this%_BASE_METHOD_ROTSYM_%finalize()
        deallocate( this%normalization_factors )

    end subroutine finalize_method_rotsym
   

    function new_wf_rotsym(m, u, coefficient) result(this)
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_ptr, c_loc
        type(_WF_ROTSYM_) :: this
        class(_METHOD_ROTSYM_), target, intent(inout) :: m
#ifdef _REAL_       
        real(kind=prec), optional, target, intent(inout) :: u(:)
        real(kind=prec), optional, intent(in) :: coefficient
#else            
        complex(kind=prec), optional, target, intent(inout) :: u(:)
        complex(kind=prec), optional, intent(in) :: coefficient
#endif
        this%_BASE_WF_ROTSYM_ = _BASE_WF_ROTSYM_(m, u, coefficient)
    end function new_wf_rotsym


end module _MODULE_
