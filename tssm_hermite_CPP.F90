#ifdef _REAL_
#define S0(x,y)  x ## _real_ ## y ## D 
#else
#define S0(x,y)  x ## _ ## y ## D 
#endif
#define S1(x,y) S0(x,y)
#define S(x) S1(x,_DIM_)
#ifdef _REAL_
 #define _WAVE_FUNCTION_ real_wave_function
 #define _COMPLEX_OR_REAL_ real
#else
 #define _WAVE_FUNCTION_ wave_function
 #define _COMPLEX_OR_REAL_ complex
#endif

#ifdef _QUADPRECISION_
#define fftw_alloc_real fftwq_alloc_real
#define fftw_alloc_complex fftwq_alloc_complex
#define fftw_free fftwq_free
#endif

#ifdef _QUADPRECISION_
module S(tssmq_hermite)
    use tssmq_base
    use tssmq_grid
    use tssmq_tensorial
    use tssmq_hermite_common
#else
module S(tssm_hermite)
    use tssm_base
    use tssm_grid
    use tssm_tensorial
    use tssm_hermite_common
#endif    
    implicit none

    private
    public :: S(hermite), S(wf_hermite)

    type, extends(S(tensorial)) :: S(hermite)
        real(kind=prec) :: gamma_x 
#if(_DIM_>=2)
        real(kind=prec) :: gamma_y 
#endif 
#if(_DIM_>=3)
        real(kind=prec) :: gamma_z 
#endif 

    contains    
        !final :: final_hermite_1D
        !! Fortran 2003 feature final seems to be not properly implemented
        !! in the gcc/gfortran compiler :(
    end type S(hermite)

    interface  S(hermite) ! constructor
        module procedure new_method
    end interface S(hermite)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    type, extends(S(wf_tensorial)) :: S(wf_hermite)
    contains
        procedure :: save
        procedure :: load
        !procedure :: evaluate
    end type S(wf_hermite)

    interface S(wf_hermite) ! constructor
        module procedure new_wf
    end interface S(wf_hermite)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains



#if(_DIM_==1)
    function new_method(nx, gamma_x ) result(this)
#elif(_DIM_==2)
    function new_method(nx, gamma_x, ny, gamma_y) result(this)
#elif(_DIM_==3)
    function new_method(nx, gamma_x, ny, gamma_y, nz, gamma_z) result(this)
#endif
        type(S(hermite)) :: this
        integer, intent(in) :: nx
        real(kind=prec), intent(in) :: gamma_x 
#if(_DIM_>=2)
        integer, intent(in) :: ny
        real(kind=prec), intent(in) :: gamma_y 
#endif
#if(_DIM_>=3)
        integer, intent(in) :: nz
        real(kind=prec), intent(in) :: gamma_z 
#endif
        integer :: i

#if(_DIM_==1)
        this%S(tensorial) = S(tensorial)(nx)
#elif(_DIM_==2)
        this%S(tensorial) = S(tensorial)(nx, ny)
#elif(_DIM_==3)
        this%S(tensorial) = S(tensorial)(nx, ny, nz)
#endif
        this%gamma_x = gamma_x
        call hermite_scaled_coeffs(nx-1, this%g%nodes_x,  this%g%weights_x, this%H1, this%gamma_x)
        this%eigenvalues1 = (/ ( -this%gamma_x*(0.5_prec + real(i, prec)), i = this%nf1min, this%nf1max) /)
#if(_DIM_>=2)
        this%gamma_y = gamma_y
        call hermite_scaled_coeffs(ny-1, this%g%nodes_y,  this%g%weights_y, this%H2, this%gamma_y)
        this%eigenvalues2 = (/ ( -this%gamma_y*(0.5_prec + real(i, prec)), i = this%nf2min, this%nf2max) /)
#endif
#if(_DIM_>=3)
        this%gamma_z = gamma_z
        call hermite_scaled_coeffs(nz-1, this%g%nodes_z,  this%g%weights_z, this%H3, this%gamma_z)
        this%eigenvalues3 = (/ ( -this%gamma_z*(0.5_prec + real(i, prec)), i = this%nf3min, this%nf3max) /)
#endif
    end function new_method


    function new_wf(m, u, coefficient) result(this)
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_ptr, c_loc
        type(S(wf_hermite)) :: this
        class(S(hermite)), target, intent(inout) :: m
#ifdef _REAL_       
#if(_DIM_==1)
        real(kind=prec), optional, target, intent(inout) :: u(:)
#elif(_DIM_==2)
        real(kind=prec), optional, target, intent(inout) :: u(:,:)
#elif(_DIM_==3)
        real(kind=prec), optional, target, intent(inout) :: u(:,:,:)
#endif
        real(kind=prec), optional, intent(in) :: coefficient
#else            
#if(_DIM_==1)
        complex(kind=prec), optional, target, intent(inout) :: u(:)
#elif(_DIM_==2)
        complex(kind=prec), optional, target, intent(inout) :: u(:,:)
#elif(_DIM_==3)
        complex(kind=prec), optional, target, intent(inout) :: u(:,:,:)
#endif
        complex(kind=prec), optional, intent(in) :: coefficient
#endif
        this%S(wf_tensorial) = S(wf_tensorial)(m, u, coefficient)
    end function new_wf


    subroutine save(this, filename, &
#ifdef _REAL_        
                    dset_name, &
#else        
                    dset_name_real, dset_name_imag, &
#endif                    
                    append)

#ifdef _NO_HDF5_
        class(S(wf_hermite)), intent(inout) :: this
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
        class(S(wf_hermite)), intent(inout) :: this
        character(len=*), intent(in) :: filename
#ifdef _REAL_        
        character(len=*), intent(in), optional :: dset_name
#else        
        character(len=*), intent(in), optional :: dset_name_real
        character(len=*), intent(in), optional :: dset_name_imag
#endif       
        logical, intent(in), optional :: append

        call this%S(wf_tensorial)%save(filename)

        select type (m=>this%m ); class is (S(hermite))
        call hdf5_write_attributes(filename, (/ "nodes_type" /), &
          (/ "hermite" /), &
          (/ character(len=0) :: /), (/ integer :: /), &
#if(_DIM_==1)        
          (/ "gamma_x" /), (/ m%gamma_x /) & 
#elif(_DIM_==2)        
          (/ "gamma_x", "gamma_y" /), (/ m%gamma_x, m%gamma_y /) & 
#elif(_DIM_==3)        
          (/ "gamma_x", "gamma_y", "gamma_z" /), (/ m%gamma_x, m%gamma_y, m%gamma_z /) & 
#endif
          )
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
        class(S(wf_hermite)), intent(inout) :: this
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
        class(S(wf_hermite)), intent(inout) :: this
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
       if (s(1:len)/="hermite") then
              stop "E: wrong nodes type"
       end if
       
       select type (m=>this%m ); class is (S(hermite))
       if ((hdf5_read_integer_attribute(filename, "nx") /= this%m%g%nx) .or. &
           (abs(hdf5_read_double_attribute(filename, "gamma_x")- m%gamma_x)>eps)) then
              stop "E: incompatible grids"
       end if       
#if(_DIM_>=2)
       if ((hdf5_read_integer_attribute(filename, "ny") /= this%m%g%ny) .or. &
           (abs(hdf5_read_double_attribute(filename, "gamma_y")- m%gamma_y)>eps)) then
              stop "E: incompatible grids"
       end if       
#endif
#if(_DIM_>=3)
       if ((hdf5_read_integer_attribute(filename, "nz") /= this%m%g%nz) .or. &
           (abs(hdf5_read_double_attribute(filename, "gamma_z")- m%gamma_z)>eps)) then
              stop "E: incompatible grids"
       end if       
#endif
       end select
    
       call this%S(wf_tensorial)%load(filename)
#endif
    end subroutine load



#ifdef _QUADPRECISION_
end module S(tssmq_hermite)
#else
end module S(tssm_hermite)
#endif
