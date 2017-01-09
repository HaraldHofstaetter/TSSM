
#ifdef _REAL_
 #ifdef _HERMITE_
   #ifdef _QUADPRECISION_
    #define _MODULE_ tssmq_generalized_laguerre_hermite_real_3d
   #else
    #define _MODULE_ tssm_generalized_laguerre_hermite_real_3d
   #endif 
    #define _METHOD_ generalized_laguerre_hermite_real_3d
    #define _WF_ wf_generalized_laguerre_hermite_real_3d
    #define _BASE_METHOD_ cylindrical_real_3d
    #define _BASE_WF_ wf_cylindrical_real_3d
 #else
   #ifdef _QUADPRECISION_
    #define _MODULE_ tssmq_generalized_laguerre_real_2d
   #else 
    #define _MODULE_ tssm_generalized_laguerre_real_2d
   #endif 
    #define _METHOD_ generalized_laguerre_real_2d
    #define _WF_ wf_generalized_laguerre_real_2d
    #define _BASE_METHOD_ polar_real_2d
    #define _BASE_WF_ wf_polar_real_2d
 #endif
 #define _COMPLEX_OR_REAL_ real
 #define _WAVE_FUNCTION_ real_wave_function
#else
 #ifdef _HERMITE_
   #ifdef _QUADPRECISION_
    #define _MODULE_ tssmq_generalized_laguerre_hermite_3d
   #else 
    #define _MODULE_ tssm_generalized_laguerre_hermite_3d
   #endif 
    #define _METHOD_ generalized_laguerre_hermite_3d
    #define _WF_ wf_generalized_laguerre_hermite_3d
    #define _BASE_METHOD_ cylindrical_3d
    #define _BASE_WF_ wf_cylindrical_3d
 #else
   #ifdef _QUADPRECISION_
    #define _MODULE_ tssmq_generalized_laguerre_2d
   #else 
    #define _MODULE_ tssm_generalized_laguerre_2d
   #endif 
    #define _METHOD_ generalized_laguerre_2d
    #define _WF_ wf_generalized_laguerre_2d
    #define _BASE_METHOD_ polar_2d
    #define _BASE_WF_ wf_polar_2d
 #endif 
 #define _COMPLEX_OR_REAL_ complex
 #define _WAVE_FUNCTION_ wave_function
#endif


module _MODULE_
#ifdef _QUADPRECISION_
    use tssmq_base
    use tssmq_grid
    use tssmq_polar
    use tssmq_fourier_common
    use tssmq_hermite_common
    use tssmq_generalized_laguerre_common
#else
    use tssm_base
    use tssm_grid
    use tssm_polar
    use tssm_fourier_common
    use tssm_hermite_common
    use tssm_generalized_laguerre_common
#endif    
    implicit none

    private
    public :: _METHOD_, _WF_

    type, extends(_BASE_METHOD_) :: _METHOD_ 
        real(kind=prec) :: gamma_r 
#ifdef _HERMITE_        
        real(kind=prec) :: gamma_z
#endif        
        real(kind=prec) :: Omega 
    contains    
        !final :: final_laguerre_1D
        !! Fortran 2003 feature final seems to be not properly implemented
        !! in the gcc/gfortran compiler :
    end type _METHOD_

    interface _METHOD_ ! constructor
        module procedure new_method
    end interface _METHOD_


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type, extends(_BASE_WF_) ::  _WF_ 
    contains
        procedure :: save
        procedure :: load
    end type _WF_ 

    interface _WF_ ! constructor
        module procedure new_wf
    end interface _WF_
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains


#ifdef _HERMITE_
    function new_method(ntheta, nfr, gamma_r, Omega, nz, gamma_z ) result(this)
#else
    function new_method(ntheta, nfr, gamma_r, Omega ) result(this)
#endif
        type(_METHOD_) :: this
        integer, intent(in) :: ntheta
        integer, intent(in) :: nfr
        real(kind=prec), intent(in) :: gamma_r 
        real(kind=prec), intent(in) :: Omega 
#ifdef _HERMITE_
        integer, intent(in) :: nz
        real(kind=prec), intent(in) :: gamma_z 
#endif
        real(kind=prec), parameter :: pi = 3.1415926535897932384626433832795028841971693993751_prec

        integer :: nr
        integer :: j, j1

        this%Omega = Omega
        this%gamma_r = gamma_r
#ifdef _HERMITE_
        this%gamma_z = gamma_z
#endif


!adjust grid parameters, TODO MPI!!!
        nr = nfr+ntheta/2

#ifndef _HERMITE_
        this%_BASE_METHOD_ = _BASE_METHOD_(ntheta, nr, nfr, .true., .true.) 
#else        
        this%_BASE_METHOD_ = _BASE_METHOD_(ntheta, nr, nfr, nz, .true., .true.) 
#endif        
          ! symmetric coefficients
          ! separated eigenvalues ...

        call generalized_laguerre_scaled_coeffs(nfr-1, ntheta, this%g%nodes_r,  this%g%weights_r, this%L, this%gamma_r )

        this%eigenvalues_r = (/ ( -this%gamma_r*(2.0_prec*real(j, prec)+1.0), j = this%nfrmin, this%nfrmax) /)

        do j = this%g%nthetamin,this%g%nthetamax
            j1 = j
            if (j1>ntheta/2) then
                j1 = j-ntheta
            end if
            this%eigenvalues_theta(j) = -this%gamma_r*abs(real(j1, kind=prec)) + this%Omega*real(j1, kind=prec)
        end do

#ifdef _HERMITE_
        call hermite_scaled_coeffs(nz, this%g%nodes_z,  this%g%weights_z, this%H_z, this%gamma_z)
        this%eigenvalues_z = (/ ( -this%gamma_z*(0.5_prec + real(j, prec)), j = this%nfzmin, this%nfzmax) /)
#endif
    end function new_method


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function new_wf(m, u, coefficient) result(this)
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_ptr, c_loc
        type(_WF_) :: this
        class(_METHOD_), target, intent(inout) :: m
#ifdef _REAL_       
#ifndef _HERMITE_
        real(kind=prec), optional, target, intent(inout) :: u(:,:)
#else
        real(kind=prec), optional, target, intent(inout) :: u(:,:,:)
#endif
        real(kind=prec), optional, intent(in) :: coefficient
#else            
#ifndef _HERMITE_
        complex(kind=prec), optional, target, intent(inout) :: u(:,:)
#else            
        complex(kind=prec), optional, target, intent(inout) :: u(:,:,:)
#endif
        complex(kind=prec), optional, intent(in) :: coefficient
#endif
        this%_BASE_WF_ = _BASE_WF_(m, u, coefficient)
    end function new_wf



    subroutine save(this, filename, &
#ifdef _REAL_        
                    dset_name, &
#else        
                    dset_name_real, dset_name_imag, &
#endif                    
                    append)

#ifdef _NO_HDF5_
        class(_WF_), intent(inout) :: this
        character(len=*), intent(in) :: filename
#ifdef _REAL_        
        character(len=*), intent(in), optional :: dset_name
#else        
        character(len=*), intent(in), optional :: dset_name_real
        character(len=*), intent(in), optional :: dset_name_imag
#endif                    
        logical, intent(in), optional :: append

        print *, "W: save not implemented"
#else
#ifdef _QUADPRECISION_
        use tssmq_hdf5
#else
        use tssm_hdf5
#endif        
        class(_WF_), intent(inout) :: this
        character(len=*), intent(in) :: filename
#ifdef _REAL_        
        character(len=*), intent(in), optional :: dset_name
#else        
        character(len=*), intent(in), optional :: dset_name_real
        character(len=*), intent(in), optional :: dset_name_imag
#endif                    
        logical, intent(in), optional :: append

        !TODO

        call this%_BASE_WF_%save(filename)

        select type (m=>this%m ); class is (_METHOD_)
#ifndef _HERMITE_
        call hdf5_write_attributes(filename, (/ "nodes_type" /), &
          (/ "generalized_laguerre" /), &
          (/ character(len=0) :: /), (/ integer :: /), &
          (/ "gamma_r" /), (/ m%gamma_r /) ) 
#else        
        call hdf5_write_attributes(filename, (/ "nodes_type" /), &
          (/ "generalized_laguerre_hermite" /), &
          (/ character(len=0) :: /), (/ integer :: /), &
          (/ "gamma_r", "gamma_z" /), (/ m%gamma_r, m%gamma_z /) ) 
#endif   
        end select
!TODO save nodes !!!
#endif        
    end subroutine save


   subroutine load(this, filename, &
#ifdef _REAL_        
                    dset_name &
#else        
                    dset_name_real, dset_name_imag &
#endif                    
                    )

#ifdef _NO_HDF5_
        class(_WF_), intent(inout) :: this
        character(len=*), intent(in) :: filename
#ifdef _REAL_        
        character(len=*), intent(in), optional :: dset_name
#else        
        character(len=*), intent(in), optional :: dset_name_real
        character(len=*), intent(in), optional :: dset_name_imag
#endif             
        print *, "W: load not implemented"
#else
#ifdef _QUADPRECISION_
        use tssmq_hdf5
#else
        use tssm_hdf5
#endif        
        class(_WF_), intent(inout) :: this
        character(len=*), intent(in) :: filename
#ifdef _REAL_        
        character(len=*), intent(in), optional :: dset_name
#else        
        character(len=*), intent(in), optional :: dset_name_real
        character(len=*), intent(in), optional :: dset_name_imag
#endif             
#ifdef _QUADPRECISION_
        real(kind=prec), parameter :: eps = epsilon(1.0_8)
#else        
        real(kind=prec), parameter :: eps = epsilon(1.0_prec)
#endif
        character(len=20) :: s
        integer :: n, len

       call hdf5_read_string_attribute(filename, "nodes_type", s, len)      
#ifndef _HERMITE_
       if (s(1:len)/="generalized_laguerre") then
              stop "E: wrong nodes type"
       end if
#else
       if (s(1:len)/="generalized_laguerre_hermite") then
              stop "E: wrong nodes type"
       end if
#endif
       select type (m=>this%m ); class is (_METHOD_)
       if (abs(hdf5_read_double_attribute(filename, "gamma_r")- m%gamma_r)>eps) then
              stop "E: incompatible grids"
       end if       
#ifdef _HERMITE_
       if (abs(hdf5_read_double_attribute(filename, "gamma_z")- m%gamma_z)>eps) then
              stop "E: incompatible grids"
       end if       
#endif
       end select

       call this%_BASE_WF_%save(filename)
#endif
    end subroutine load


end module _MODULE_ 
