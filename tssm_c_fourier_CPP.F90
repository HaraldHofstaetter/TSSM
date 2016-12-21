#ifdef _REAL_
#define S0(x,y)  x ## _real_ ## y ## d 
#else
#define S0(x,y)  x ## _ ## y ## d 
#endif
#define S1(x,y) S0(x,y)
#define S(x) S1(x,_DIM_)
!#define SC0(x) #x
!#define SC1(x) SC0(x)
!#define SC(x) SC1(S(x))

#ifdef _QUADPRECISION_
    #define C0(x) tssmq_ ## x 
#else    
    #define C0(x) tssm_ ## x
#endif    
#define C1(x) C0(x)
#define C(x) C1(S(x))
#define SC0(x) #x
#define SC1(x) SC0(x)
#define SC(x) SC1(C(x))


#ifdef _REAL_
 #define _WAVE_FUNCTION_ real_wave_function
 #define _COMPLEX_OR_REAL_ real
#else
 #define _WAVE_FUNCTION_ wave_function
 #define _COMPLEX_OR_REAL_ complex
#endif



#ifdef _QUADPRECISION_
module S(tssmq_c_fourier)
    use tssmq_base
    use tssmq_fourier
#else
module S(tssm_c_fourier)
    use tssm_base
    use tssm_fourier
#endif
    implicit none
    private

contains

!    subroutine c_initialize_tssm_fourier() bind(c)
!        call initialize_tssm_fourier
!    end subroutine c_initialize_tssm_fourier 

    function c_new(nx, xmin, xmax, &
#if(_DIM_>=2)
                              ny, ymin, ymax, &
#endif
#if(_DIM_>=3)
                              nz, zmin, zmax, &
#endif
                              boundary_conditions) &
        result(this) bind(c, name=SC(new_fourier))
        use iso_c_binding
        type(c_ptr) :: this 
        real(kind=prec),  value :: xmin 
        real(kind=prec),  value :: xmax
        integer, value :: nx
#if(_DIM_>=2)
        real(kind=prec),  value :: ymin 
        real(kind=prec),  value :: ymax
        integer, value :: ny
#endif
#if(_DIM_>=3)
        real(kind=prec),  value :: zmin 
        real(kind=prec),  value :: zmax
        integer, value :: nz
#endif
        integer, value :: boundary_conditions

        type(S(fourier)), pointer :: m
        allocate( m )
        m = S(fourier)(nx, xmin, xmax, &
#if(_DIM_>=2)
                       ny, ymin, ymax, &
#endif
#if(_DIM_>=3)
                       nz, zmin, zmax, &
#endif
                       boundary_conditions)
        this =  c_loc(m)
    end function c_new


    subroutine c_finalize(m) &
        bind(c, name=SC(finalize_fourier))
        use iso_c_binding
        type(c_ptr), value :: m
        type(S(fourier)), pointer :: mp

        call c_f_pointer(m, mp)
        call mp%finalize
        deallocate( mp )
    end subroutine c_finalize

    subroutine c_set_propagate_time_together_with_A(m, flag) &
        bind(c, name=SC(set_propagate_time_together_with_A_fourier)) 
        use iso_c_binding
        type(c_ptr), value :: m 
        logical, value :: flag 
        type(S(fourier)), pointer :: mp

        call c_f_pointer(m, mp)
        mp%propagate_time_together_with_A = flag 
    end subroutine c_set_propagate_time_together_with_A

    function c_get_propagate_time_together_with_A(m) &
        result(flag) bind(c, name=SC(get_propagate_time_together_with_A_fourier))
        use iso_c_binding
        type(c_ptr), value :: m 
        type(S(fourier)), pointer :: mp
        logical :: flag 

        call c_f_pointer(m, mp)
        flag = mp%propagate_time_together_with_A
    end function c_get_propagate_time_together_with_A

    function c_new_wf(m) &
        result(this) bind(c, name=SC(new_wf_fourier))
        use iso_c_binding
        type(c_ptr) :: this 
        type(c_ptr), value :: m 
        type(S(fourier)), pointer :: mp
        type(S(wf_fourier)), pointer :: psi

        call c_f_pointer(m, mp)
        allocate( psi )
        psi = S(wf_fourier)(mp)
        this =  c_loc(psi)
    end function c_new_wf

    subroutine c_finalize_wf(psi) &
        bind(c, name=SC(finalize_wf_fourier))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%finalize
        deallocate( psip )
    end subroutine c_finalize_wf


    function c_is_real_space_wf(psi) &
        result(ans) bind(c, name=SC(is_real_space_wf_fourier))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip
        logical :: ans

        call c_f_pointer(psi, psip)
        ans = psip%is_real_space()
    end function c_is_real_space_wf

    subroutine c_to_real_space_wf(psi) &
        bind(c, name=SC(to_real_space_wf_fourier)) 
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%to_real_space
    end subroutine c_to_real_space_wf
  
    subroutine c_to_frequency_space_wf(psi) &
        bind(c, name=SC(to_frequency_space_wf_fourier)) 
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%to_frequency_space
    end subroutine c_to_frequency_space_wf

    subroutine c_set_time_wf(psi, t) &
        bind(c, name=SC(set_time_wf_fourier)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: t
        type(S(wf_fourier)), pointer :: psip

        call c_f_pointer(psi, psip)
        psip%time = t
    end subroutine c_set_time_wf

    function c_get_time_wf(psi) &
        result(ans) bind(c, name=SC(get_time_wf_fourier))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip
        _COMPLEX_OR_REAL_(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        ans = psip%time
    end function c_get_time_wf

    subroutine c_propagate_time_wf(psi, dt) &
        bind(c, name=SC(propagate_time_wf_fourier)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(S(wf_fourier)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%propagate_time(dt)
    end subroutine c_propagate_time_wf
    

    subroutine c_propagate_A_wf(psi, dt) &
        bind(c, name=SC(propagate_A_wf_fourier)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(S(wf_fourier)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%propagate_A(dt)
    end subroutine c_propagate_A_wf


    subroutine c_propagate_A_derivative_wf(psi, wf, dt) & 
        bind(c, name=SC(propagate_A_derivative_wf_fourier)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        type(c_ptr), value :: wf
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(S(wf_fourier)), pointer :: psip
        type(S(wf_fourier)), pointer :: wfp 

        call c_f_pointer(psi, psip)
        call c_f_pointer(wf, wfp)
        call psip%propagate_A_derivative(wfp, dt)
    end subroutine c_propagate_A_derivative_wf
    

    subroutine c_propagate_B_wf(psi, dt) & 
        bind(c, name=SC(propagate_B_wf_fourier)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(S(wf_fourier)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%propagate_B(dt)
    end subroutine c_propagate_B_wf

    subroutine c_propagate_C_wf(psi, dt) &
        bind(c, name=SC(propagate_C_wf_fourier)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(S(wf_fourier)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%propagate_C(dt)
    end subroutine c_propagate_C_wf

    subroutine c_add_apply_A_wf(this, other, coefficient) &
        bind(c, name=SC(add_apply_A_wf_fourier)) 
        use iso_c_binding
        type(c_ptr), value :: this
        type(c_ptr), value :: other 
        _COMPLEX_OR_REAL_(kind=prec), value :: coefficient
        type(S(wf_fourier)), pointer :: thisp
        type(S(wf_fourier)), pointer :: otherp

        call c_f_pointer(this, thisp)
        call c_f_pointer(other, otherp)
        call thisp%add_apply_A(otherp, coefficient)
    end subroutine c_add_apply_A_wf


    subroutine c_save_wf(psi, filename, filename_length) &
        bind(c, name=SC(save_wf_fourier)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        integer(c_int), intent(in), value :: filename_length
        character(kind=c_char), dimension(filename_length), intent(in):: filename
        character(len=filename_length, kind=c_char) :: fn 
        type(S(wf_fourier)), pointer :: psip
        integer :: j
        call c_f_pointer(psi, psip)
        do j=1,filename_length
            fn(j:j) = filename(j)
        end do
        call psip%save(fn)
    end subroutine c_save_wf

    subroutine c_load_wf(psi, filename, filename_length) &
        bind(c, name=SC(load_wf_fourier)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        integer(c_int), intent(in), value :: filename_length
        character(kind=c_char), dimension(filename_length), intent(in):: filename
        character(len=filename_length, kind=c_char) :: fn 
        type(S(wf_fourier)), pointer :: psip
        integer :: j
        call c_f_pointer(psi, psip)
        do j=1,filename_length
            fn(j:j) = filename(j)
        end do
        call psip%load(fn)
    end subroutine c_load_wf

    function c_norm_wf(psi) &
        result(ans) bind(c, name=SC(norm_wf_fourier))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip
        real(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        ans = psip%norm()
    end function c_norm_wf

    function c_norm_in_frequency_space_wf(psi) &
        result(ans) bind(c, name=SC(norm_in_frequency_space_wf_fourier))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip
        real(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        ans = psip%norm_in_frequency_space()
    end function c_norm_in_frequency_space_wf


    function c_inner_product_wf(psi, other) &
        result(ans) bind(c, name=SC(inner_product_wf_fourier))
        use iso_c_binding
        type(c_ptr), value :: psi
        type(c_ptr), value :: other
        type(S(wf_fourier)), pointer :: psip
        type(S(wf_fourier)), pointer :: otherp
        _COMPLEX_OR_REAL_(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        call c_f_pointer(other, otherp)
        ans = psip%inner_product(otherp)
   end function c_inner_product_wf
    

   function c_normalize_wf(psi) &
        result(ans) bind(c, name=SC(normalize_wf_fourier))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip
        real(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        call psip%normalize(ans)
    end function c_normalize_wf

   subroutine c_scale_wf(psi, factor) &
        bind(c, name=SC(scale_wf_fourier))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip
        _COMPLEX_OR_REAL_(kind=prec), value :: factor 

        call c_f_pointer(psi, psip)
        call psip%scale(factor)
    end subroutine c_scale_wf

   subroutine c_axpy_wf(this, other, factor) &
        bind(c, name=SC(axpy_wf_fourier)) 
        use iso_c_binding
        type(c_ptr), value :: this
        type(c_ptr), value :: other 
        _COMPLEX_OR_REAL_(kind=prec), value :: factor
        type(S(wf_fourier)), pointer :: thisp
        type(S(wf_fourier)), pointer :: otherp

        call c_f_pointer(this, thisp)
        call c_f_pointer(other, otherp)
        call thisp%axpy(otherp, factor)
    end subroutine c_axpy_wf



    function c_get_data_wf(psi, dim) &
        result(up) bind(c, name=SC(get_data_wf_fourier))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(c_ptr) :: up
        integer, intent(out)  :: dim(_DIM_)
        type(S(wf_fourier)), pointer :: psip

        call c_f_pointer(psi, psip)
        !up = c_loc(psip%u)
!#if(_DIM_==1)
!        up = c_loc(psip%u(psip%m%g%m1min))
!#elif(_DIM_==2)
!        up = c_loc(psip%u(psip%m%g%m1min, psip%m%g%m2min))
!#elif(_DIM_==3)
!        up = c_loc(psip%u(psip%m%g%m1min, psip%m%g%m2min, psip%m%g%m3min))
!#endif
        up = c_loc(psip%up(1))
        dim(1) = psip%m%g%m1max-psip%m%g%m1min+1
#if(_DIM_>=2)        
        dim(2) = psip%m%g%m2max-psip%m%g%m2min+1
#endif        
#if(_DIM_>=3)        
        dim(3) = psip%m%g%m3max-psip%m%g%m3min+1
#endif        
    end function c_get_data_wf


    function c_get_eigenvalues(m, dim, which) &
        result(evp) bind(c, name=SC(get_eigenvalues_fourier))
        use iso_c_binding
        type(c_ptr), value :: m
        type(c_ptr) :: evp
        integer, intent(out)  :: dim(1)
        integer, value :: which
        type(S(fourier)), pointer :: mp

        call c_f_pointer(m, mp)
        select case (which)
        case (1)
           !evp = c_loc(mp%eigenvalues1
           evp = c_loc(mp%eigenvalues1(mp%nf1min))
           dim(1) = mp%nf1max-mp%nf1min+1
#if(_DIM_>=2)        
        case (2)
           !evp = c_loc(mp%eigenvalues2)
           evp = c_loc(mp%eigenvalues2(mp%nf2min))
           dim(1) = mp%nf2max-mp%nf2min+1
#endif        
#if(_DIM_>=3)        
        case (3)
           !evp = c_loc(mp%eigenvalues3)
           evp = c_loc(mp%eigenvalues3(mp%nf3min))
           dim(1) = mp%nf3max-mp%nf3min+1
#endif        
        end select
    end function c_get_eigenvalues


    function c_get_nodes(m, dim, which) &
        result(np) bind(c, name=SC(get_nodes_fourier))
        use iso_c_binding
        type(c_ptr), value :: m
        type(c_ptr) :: np
        integer, intent(out)  :: dim(1)
        integer, value :: which
        type(S(fourier)), pointer :: mp

        call c_f_pointer(m, mp)
        select case (which)
        case (1)
            np = c_loc(mp%g%nodes_x)
            dim(1) = mp%g%n1max-mp%g%n1min+1
#if(_DIM_>=2)        
            np = c_loc(mp%g%nodes_y)
            dim(1) = mp%g%n2max-mp%g%n2min+1
#endif
#if(_DIM_>=3)        
            np = c_loc(mp%g%nodes_z)
            dim(1) = mp%g%n3max-mp%g%n3min+1
#endif
        end select
    end function c_get_nodes
    
    function c_get_nx(m) &
        result(nx) bind(c, name=SC(get_nx_fourier))
        use iso_c_binding
        type(c_ptr), value :: m
        integer :: nx
        type(S(fourier)), pointer :: mp
        call c_f_pointer(m, mp)
        nx = mp%g%nx
    end function c_get_nx   


    function c_get_xmin(m) &
        result(xmin) bind(c, name=SC(get_xmin_fourier))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: xmin
        type(S(fourier)), pointer :: mp
        call c_f_pointer(m, mp)
        xmin = mp%g%xmin
    end function c_get_xmin   

    function c_get_xmax(m) &
        result(xmax) bind(c, name=SC(get_xmax_fourier))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: xmax
        type(S(fourier)), pointer :: mp
        call c_f_pointer(m, mp)
        xmax = mp%g%xmax
    end function c_get_xmax     

#if(_DIM_>=2)
    function c_get_ny(m) &
        result(ny) bind(c, name=SC(get_ny_fourier))
        use iso_c_binding
        type(c_ptr), value :: m
        integer :: ny
        type(S(fourier)), pointer :: mp
        call c_f_pointer(m, mp)
        ny = mp%g%ny
    end function c_get_ny   


    function c_get_ymin(m) &
        result(ymin) bind(c, name=SC(get_ymin_fourier))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: ymin
        type(S(fourier)), pointer :: mp
        call c_f_pointer(m, mp)
        ymin = mp%g%ymin
    end function c_get_ymin   

    function c_get_ymax(m) &
        result(ymax) bind(c, name=SC(get_ymax_fourier))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: ymax
        type(S(fourier)), pointer :: mp
        call c_f_pointer(m, mp)
        ymax = mp%g%ymax
    end function c_get_ymax     

#endif

#if(_DIM_>=3)
    function c_get_nz(m) &
        result(nz) bind(c, name=SC(get_nz_fourier))
        use iso_c_binding
        type(c_ptr), value :: m
        integer :: nz
        type(S(fourier)), pointer :: mp
        call c_f_pointer(m, mp)
        nz = mp%g%nz
    end function c_get_nz   


    function c_get_zmin(m) &
        result(zmin) bind(c, name=SC(get_zmin_fourier))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: zmin
        type(S(fourier)), pointer :: mp
        call c_f_pointer(m, mp)
        zmin = mp%g%zmin
    end function c_get_zmin   

    function c_get_zmax(m) &
        result(zmax) bind(c, name=SC(get_zmax_fourier))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: zmax
        type(S(fourier)), pointer :: mp
        call c_f_pointer(m, mp)
        zmax = mp%g%zmax
    end function c_get_zmax     
#endif
    

#ifndef _REAL_ 
    subroutine c_rset_wf(psi, f) &
        bind(c, name=SC(rset_wf_fourier))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip
        interface 
#if(_DIM_==1)
           function f(x) bind(c)
#elif(_DIM_==2)
           function f(x, y) bind(c)
#elif(_DIM_==3)
           function f(x, y, z) bind(c)
#endif           
               import prec
               real(kind=prec), value :: x
#if(_DIM_>=2)
               real(kind=prec), value :: y
#endif               
#if(_DIM_>=3)
               real(kind=prec), value :: z
#endif               
               real(kind=prec) :: f 
           end function f
        end interface 

        call c_f_pointer(psi, psip)
        call psip%rset(ff)
    contains
#if(_DIM_==1)
        function ff(x) 
#elif(_DIM_==2)
        function ff(x, y)
#elif(_DIM_==3)
        function ff(x, y, z)
#endif          
            real(kind=prec), intent(in) :: x
#if(_DIM_>=2)
            real(kind=prec), intent(in) :: y
#endif               
#if(_DIM_>=3)
            real(kind=prec), intent(in) :: z
#endif               
            real(kind=prec) :: ff 
#if(_DIM_==1)
            ff = f(x)
#elif(_DIM_==2)
            ff = f(x, y)
#elif(_DIM_==3)
            ff = f(x, y, z)
#endif
        end function ff
    end subroutine c_rset_wf
#endif    

    subroutine c_set_wf(psi, f) &
        bind(c, name=SC(set_wf_fourier))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip
        interface 
#if(_DIM_==1)
           function f(x) bind(c)
#elif(_DIM_==2)
           function f(x, y) bind(c)
#elif(_DIM_==3)
           function f(x, y, z) bind(c)
#endif           
               import prec
               real(kind=prec), value :: x
#if(_DIM_>=2)
               real(kind=prec), value :: y
#endif               
#if(_DIM_>=3)
               real(kind=prec), value :: z
#endif               
               _COMPLEX_OR_REAL_(kind=prec) :: f 
           end function f
        end interface 

        call c_f_pointer(psi, psip)
        call psip%set(ff)
    contains
#if(_DIM_==1)
        function ff(x) 
#elif(_DIM_==2)
        function ff(x, y)
#elif(_DIM_==3)
        function ff(x, y, z)
#endif          
            real(kind=prec), intent(in) :: x
#if(_DIM_>=2)
            real(kind=prec), intent(in) :: y
#endif               
#if(_DIM_>=3)
            real(kind=prec), intent(in) :: z
#endif               
            _COMPLEX_OR_REAL_(kind=prec) :: ff 
#if(_DIM_==1)
            ff = f(x)
#elif(_DIM_==2)
            ff = f(x, y)
#elif(_DIM_==3)
            ff = f(x, y, z)
#endif
        end function ff
    end subroutine c_set_wf


#ifndef _REAL_ 
    subroutine c_rset_t_wf(psi, f, t) &
        bind(c, name=SC(rset_t_wf_fourier))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip
        real(kind=prec), value :: t
        interface 
#if(_DIM_==1)
           function f(x, t) bind(c)
#elif(_DIM_==2)
           function f(x, y, t) bind(c)
#elif(_DIM_==3)
           function f(x, y, z, t) bind(c)
#endif           
               import prec
               real(kind=prec), value :: x
#if(_DIM_>=2)
               real(kind=prec), value :: y
#endif               
#if(_DIM_>=3)
               real(kind=prec), value :: z
#endif               
               real(kind=prec) :: f 
               real(kind=prec), value :: t
           end function f
        end interface 

        call c_f_pointer(psi, psip)
        call psip%rset_t(ff, t)
    contains
#if(_DIM_==1)
        function ff(x, t) 
#elif(_DIM_==2)
        function ff(x, y, t)
#elif(_DIM_==3)
        function ff(x, y, z, t)
#endif          
            real(kind=prec), intent(in) :: x
#if(_DIM_>=2)
            real(kind=prec), intent(in) :: y
#endif               
#if(_DIM_>=3)
            real(kind=prec), intent(in) :: z
#endif               
            real(kind=prec), intent(in) :: t
            real(kind=prec) :: ff 
#if(_DIM_==1)
            ff = f(x, t)
#elif(_DIM_==2)
            ff = f(x, y, t)
#elif(_DIM_==3)
            ff = f(x, y, z, t)
#endif
        end function ff
    end subroutine c_rset_t_wf
#endif    

    subroutine c_set_t_wf(psi, f, t) &
        bind(c, name=SC(set_t_wf_fourier))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip
        real(kind=prec), value :: t
        interface 
#if(_DIM_==1)
           function f(x, t) bind(c)
#elif(_DIM_==2)
           function f(x, y, t) bind(c)
#elif(_DIM_==3)
           function f(x, y, z, t) bind(c)
#endif           
               import prec
               real(kind=prec), value :: x
#if(_DIM_>=2)
               real(kind=prec), value :: y
#endif               
#if(_DIM_>=3)
               real(kind=prec), value :: z
#endif               
               _COMPLEX_OR_REAL_(kind=prec) :: f 
               real(kind=prec), value :: t
           end function f
        end interface 

        call c_f_pointer(psi, psip)
        call psip%set_t(ff, t)
    contains
#if(_DIM_==1)
        function ff(x, t) 
#elif(_DIM_==2)
        function ff(x, y, t)
#elif(_DIM_==3)
        function ff(x, y, z, t)
#endif          
            real(kind=prec), intent(in) :: x
#if(_DIM_>=2)
            real(kind=prec), intent(in) :: y
#endif               
#if(_DIM_>=3)
            real(kind=prec), intent(in) :: z
#endif               
            real(kind=prec), intent(in) :: t
            _COMPLEX_OR_REAL_(kind=prec) :: ff 
#if(_DIM_==1)
            ff = f(x, t)
#elif(_DIM_==2)
            ff = f(x, y, t)
#elif(_DIM_==3)
            ff = f(x, y, z, t)
#endif
        end function ff
    end subroutine c_set_t_wf



    subroutine c_copy_wf(psi, source) &
        bind(c, name=SC(copy_wf_fourier))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(c_ptr), value :: source 
        type(S(wf_fourier)), pointer :: psip, sourcep

        call c_f_pointer(psi, psip)
        call c_f_pointer(source, sourcep)
        call psip%copy(sourcep)
    end subroutine c_copy_wf

    function c_distance_wf(psi1, psi2) &
        result(d) bind(c, name=SC(distance_wf_fourier))
        use iso_c_binding
        type(c_ptr), value :: psi1 
        type(c_ptr), value :: psi2 
        real(kind = prec) :: d
        type(S(wf_fourier)), pointer :: psi1p, psi2p

        call c_f_pointer(psi1, psi1p)
        call c_f_pointer(psi2, psi2p)
        d =  psi1p%distance(psi2p)
    end function c_distance_wf


#ifdef _QUADPRECISION_
end module S(tssmq_c_fourier)
#else
end module S(tssm_c_fourier)
#endif
