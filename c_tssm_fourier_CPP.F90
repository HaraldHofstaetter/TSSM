#ifdef _REAL_
#define S0(x,y)  x ## _real_ ## y ## D 
#define W0(x,y)  x ## _real_ ## y ## D_w 
#else
#define S0(x,y)  x ## _ ## y ## D 
#define W0(x,y)  x ## _ ## y ## D_w
#endif
#define S1(x,y) S0(x,y)
#define S(x) S1(x,_DIM_)
#ifdef _QUADPRECISION_
#define W1(x,y) W0(x,y)
#define W(x) W1(x,_DIM_)
#else
#define W(x) S1(x,_DIM_)
#endif
#ifdef _REAL_
 #define _WAVE_FUNCTION_ real_wave_function
 #define _COMPLEX_OR_REAL_ real
#else
 #define _WAVE_FUNCTION_ wave_function
 #define _COMPLEX_OR_REAL_ complex
#endif



module S(c_tssm_fourier)
    use tssm_fourier
    implicit none

contains

!    subroutine c_initialize_tssm_fourier() bind(c)
!        call initialize_tssm_fourier
!    end subroutine c_initialize_tssm_fourier 

    function W(c_new_fourier)(nx, xmin, xmax, &
#if(_DIM_>=2)
                              ny, ymin, ymax, &
#endif
#if(_DIM_>=3)
                              nz, zmin, zmax, &
#endif
                              boundary_conditions) &
        result(this) bind(c)
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
    end function W(c_new_fourier)


    subroutine W(c_finalize_fourier)(m)  bind(c)
        use iso_c_binding
        type(c_ptr) :: m
        type(S(fourier)), pointer :: mp

        call c_f_pointer(m, mp)
        call mp%finalize
        deallocate( mp )
    end subroutine W(c_finalize_fourier)



    function W(c_new_wf_fourier)(m) &
        result(this) bind(c)
        use iso_c_binding
        type(c_ptr) :: this 
        type(c_ptr), value :: m 
        type(S(fourier)), pointer :: mp
        type(S(wf_fourier)), pointer :: psi

        call c_f_pointer(m, mp)
        allocate( psi )
        psi = S(wf_fourier)(mp)
        this =  c_loc(psi)
    end function W(c_new_wf_fourier)

    subroutine W(c_finalize_wf_fourier)(psi)  bind(c)
        use iso_c_binding
        type(c_ptr) :: psi 
        type(S(wf_fourier)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%finalize
        deallocate( psip )
    end subroutine W(c_finalize_wf_fourier)


    function W(c_is_real_space_wf_fourier)(psi) &
        result(ans) bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip
        logical :: ans

        call c_f_pointer(psi, psip)
        ans = psip%is_real_space
    end function W(c_is_real_space_wf_fourier)

    subroutine W(c_to_real_space_wf_fourier)(psi) bind(c) 
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%to_real_space
    end subroutine W(c_to_real_space_wf_fourier)
  
    subroutine W(c_to_frequency_space_wf_fourier)(psi) bind(c) 
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%to_frequency_space
    end subroutine W(c_to_frequency_space_wf_fourier)

    subroutine W(c_propagate_A_wf_fourier)(psi, dt) bind(c) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(S(wf_fourier)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%propagate_A(dt)
    end subroutine W(c_propagate_A_wf_fourier)

    subroutine W(c_propagate_B_wf_fourier)(psi, dt) bind(c) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(S(wf_fourier)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%propagate_B(dt)
    end subroutine W(c_propagate_B_wf_fourier)

    subroutine W(c_propagate_C_wf_fourier)(psi, dt) bind(c) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(S(wf_fourier)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%propagate_C(dt)
    end subroutine W(c_propagate_C_wf_fourier)

    subroutine W(c_add_apply_A_wf_fourier)(this, other, coefficient) bind(c) 
        use iso_c_binding
        type(c_ptr), value :: this
        type(c_ptr), value :: other 
        _COMPLEX_OR_REAL_(kind=prec), value :: coefficient
        type(S(wf_fourier)), pointer :: thisp
        type(S(wf_fourier)), pointer :: otherp

        call c_f_pointer(this, thisp)
        call c_f_pointer(other, otherp)
        call thisp%add_apply_A(otherp, coefficient)
    end subroutine W(c_add_apply_A_wf_fourier)


    subroutine W(c_save_wf_fourier)(psi, filename, filename_length) bind(c) 
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
    end subroutine W(c_save_wf_fourier)

    subroutine W(c_load_wf_fourier)(psi, filename, filename_length) bind(c) 
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
    end subroutine W(c_load_wf_fourier)

    function W(c_norm2_wf_fourier)(psi) &
        result(ans) bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip
        real(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        ans = psip%norm2()
    end function W(c_norm2_wf_fourier)

    function W(c_norm2_in_frequency_space_wf_fourier)(psi) &
        result(ans) bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip
        real(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        ans = psip%norm2_in_frequency_space()
    end function W(c_norm2_in_frequency_space_wf_fourier)
    

   function W(c_normalize_wf_fourier)(psi) &
        result(ans) bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip
        real(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        call psip%normalize(ans)
    end function W(c_normalize_wf_fourier)

   subroutine W(c_scale_wf_fourier)(psi, factor) bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_fourier)), pointer :: psip
        _COMPLEX_OR_REAL_(kind=prec), value :: factor 

        call c_f_pointer(psi, psip)
        call psip%scale(factor)
    end subroutine W(c_scale_wf_fourier)

   subroutine W(c_axpy_wf_fourier)(this, other, factor) bind(c) 
        use iso_c_binding
        type(c_ptr), value :: this
        type(c_ptr), value :: other 
        _COMPLEX_OR_REAL_(kind=prec), value :: factor
        type(S(wf_fourier)), pointer :: thisp
        type(S(wf_fourier)), pointer :: otherp

        call c_f_pointer(this, thisp)
        call c_f_pointer(other, otherp)
        call thisp%axpy(otherp, factor)
    end subroutine W(c_axpy_wf_fourier)



    function W(c_get_data_wf_fourier)(psi, dim) &
        result(up) bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(c_ptr)  :: up
        integer, intent(out)  :: dim(_DIM_)
        type(S(wf_fourier)), pointer :: psip

        call c_f_pointer(psi, psip)
        up = c_loc(psip%u)
        dim(1) = psip%m%g%m1max-psip%m%g%m1min+1
#if(_DIM_>=2)        
        dim(2) = psip%m%g%m2max-psip%m%g%m2min+1
#endif        
#if(_DIM_>=3)        
        dim(3) = psip%m%g%m3max-psip%m%g%m3min+1
#endif        
    end function W(c_get_data_wf_fourier)


    function W(c_get_eigenvalues_fourier)(m, dim, which) &
        result(evp) bind(c)
        use iso_c_binding
        type(c_ptr), value :: m
        type(c_ptr) :: evp
        integer, intent(out)  :: dim(1)
        integer, value :: which
        type(S(fourier)), pointer :: mp

        call c_f_pointer(m, mp)
        select case (which)
        case (1)
           evp = c_loc(mp%eigenvalues1)
           dim(1) = mp%nf1max-mp%nf1min+1
#if(_DIM_>=2)        
        case (2)
           evp = c_loc(mp%eigenvalues2)
           dim(1) = mp%nf2max-mp%nf2min+1
#endif        
#if(_DIM_>=3)        
        case (3)
           evp = c_loc(mp%eigenvalues3)
           dim(1) = mp%nf3max-mp%nf3min+1
#endif        
        end select
    end function W(c_get_eigenvalues_fourier)


    function W(c_get_nodes_fourier)(m, dim, which) &
        result(np) bind(c)
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
    end function W(c_get_nodes_fourier)
    
    function W(c_get_nx_fourier)(m) &
        result(nx) bind(c)
        use iso_c_binding
        type(c_ptr), value :: m
        integer :: nx
        type(S(fourier)), pointer :: mp
        call c_f_pointer(m, mp)
        nx = mp%g%nx
    end function W(c_get_nx_fourier)   


    function W(c_get_xmin_fourier)(m) &
        result(xmin) bind(c)
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: xmin
        type(S(fourier)), pointer :: mp
        call c_f_pointer(m, mp)
        xmin = mp%g%xmin
    end function W(c_get_xmin_fourier)   

    function W(c_get_xmax_fourier)(m) &
        result(xmax) bind(c)
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: xmax
        type(S(fourier)), pointer :: mp
        call c_f_pointer(m, mp)
        xmax = mp%g%xmax
    end function W(c_get_xmax_fourier)     

#if(_DIM_>=2)
    function W(c_get_ny_fourier)(m) &
        result(ny) bind(c)
        use iso_c_binding
        type(c_ptr), value :: m
        integer :: ny
        type(S(fourier)), pointer :: mp
        call c_f_pointer(m, mp)
        ny = mp%g%ny
    end function W(c_get_ny_fourier)   


    function W(c_get_ymin_fourier)(m) &
        result(ymin) bind(c)
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: ymin
        type(S(fourier)), pointer :: mp
        call c_f_pointer(m, mp)
        ymin = mp%g%ymin
    end function W(c_get_ymin_fourier)   

    function W(c_get_ymax_fourier)(m) &
        result(ymax) bind(c)
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: ymax
        type(S(fourier)), pointer :: mp
        call c_f_pointer(m, mp)
        ymax = mp%g%ymax
    end function W(c_get_ymax_fourier)     

#endif

#if(_DIM_>=3)
    function W(c_get_nz_fourier)(m) &
        result(nz) bind(c)
        use iso_c_binding
        type(c_ptr), value :: m
        integer :: nz
        type(S(fourier)), pointer :: mp
        call c_f_pointer(m, mp)
        nz = mp%g%nz
    end function W(c_get_nz_fourier)   


    function W(c_get_zmin_fourier)(m) &
        result(zmin) bind(c)
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: zmin
        type(S(fourier)), pointer :: mp
        call c_f_pointer(m, mp)
        zmin = mp%g%zmin
    end function W(c_get_zmin_fourier)   

    function W(c_get_zmax_fourier)(m) &
        result(zmax) bind(c)
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: zmax
        type(S(fourier)), pointer :: mp
        call c_f_pointer(m, mp)
        zmax = mp%g%zmax
    end function W(c_get_zmax_fourier)     
#endif
    

#ifndef _REAL_ 
    subroutine W(c_rset_wf_fourier)(psi, f) &
        bind(c)
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
    end subroutine W(c_rset_wf_fourier)
#endif    

    subroutine W(c_set_wf_fourier)(psi, f) &
        bind(c)
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
    end subroutine W(c_set_wf_fourier)


#ifndef _REAL_ 
    subroutine W(c_rset_t_wf_fourier)(psi, f, t) &
        bind(c)
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
    end subroutine W(c_rset_t_wf_fourier)
#endif    

    subroutine W(c_set_t_wf_fourier)(psi, f, t) &
        bind(c)
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
    end subroutine W(c_set_t_wf_fourier)



    subroutine W(c_copy_wf_fourier)(psi, source) &
        bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(c_ptr), value :: source 
        type(S(wf_fourier)), pointer :: psip, sourcep

        call c_f_pointer(psi, psip)
        call c_f_pointer(source, sourcep)
        call psip%copy(sourcep)
    end subroutine W(c_copy_wf_fourier)

    function W(c_distance_wf_fourier)(psi1, psi2) &
        result(d) bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi1 
        type(c_ptr), value :: psi2 
        real(kind = prec) :: d
        type(S(wf_fourier)), pointer :: psi1p, psi2p

        call c_f_pointer(psi1, psi1p)
        call c_f_pointer(psi2, psi2p)
        d =  psi1p%distance(psi2p)
    end function W(c_distance_wf_fourier)


end module S(c_tssm_fourier)
