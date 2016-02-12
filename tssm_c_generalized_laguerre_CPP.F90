#if 0
!!! TODO TODO TODO !!!


#ifdef _HERMITE_
    #define _DIM_ 3
#else    
    #define _DIM_ 2
#endif 

#ifdef _REAL_
 #ifdef _HERMITE_
    #define _MODULE_ tssm_generalized_laguerre_hermite_real_3d
    #define _C_MODULE_ c_tssm_generalized_laguerre_hermite_real_3d
    #define _METHOD_ generalized_laguerre_hermite_real_3d
    #define _WF_ wf_generalized_laguerre_hermite_real_3d
    #define S(x) x ## _gen_laguerre_hermite_real_3d
 #else
    #define _MODULE_ tssm_generalized_laguerre_real_2d
    #define _C_MODULE_ c_tssm_generalized_laguerre_real_2d
    #define _METHOD_ generalized_laguerre_real_2d
    #define _WF_ wf_generalized_laguerre_real_2d
    #define S(x) x ## _gen_laguerre_real_2d
 #endif
 #define _WAVE_FUNCTION_ real_wave_function
 #define _COMPLEX_OR_REAL_ real
#else
 #ifdef _HERMITE_
    #define _MODULE_ tssm_generalized_laguerre_hermite_3d
    #define _C_MODULE_ c_tssm_generalized_laguerre_hermite_3d
    #define _METHOD_ generalized_laguerre_hermite_3d
    #define _WF_ wf_generalized_laguerre_hermite_3d
    #define S(x) x ## _gen_laguerre_hermite_3d
 #else
    #define _MODULE_ tssm_generalized_laguerre_2d
    #define _C_MODULE_ c_tssm_generalized_laguerre_2d
    #define _METHOD_ generalized_laguerre_2d
    #define _WF_ wf_generalized_laguerre_2d
    #define S(x) x ## _gen_laguerre_2d
 #endif 
 #define _WAVE_FUNCTION_ wave_function
 #define _COMPLEX_OR_REAL_ complex
#endif

module _C_MODULE_ 
    use tssm_generalized_laguerre
    implicit none

contains

#ifdef _HERMITE_
    function S(c_new)(M, K, gamma_r, Omega, nz, gamma_z) &
#else
    function S(c_new)(M, K, gamma_r, Omega) &
#endif
        result(this) bind(c)
        use iso_c_binding
        type(c_ptr) :: this 
        integer, value :: M
        integer, value :: K
        real(kind=prec),  value :: gamma_r 
        real(kind=prec),  value :: Omega 
#ifdef _HERMITE_
        integer, value :: nz 
        real(kind=prec),  value :: gamma_z 
#endif

        type(_METHOD_), pointer :: meth
        allocate( meth )
#ifdef _HERMITE_
        meth = _METHOD_( M, K, gamma_r, Omega, nz, gamma_z )
#else
        meth = _METHOD_( M, K, gamma_r, Omega )
#endif
        this =  c_loc(meth)
    end function S(c_new)

    subroutine S(c_finalize)(m)  bind(c)
        use iso_c_binding
        type(c_ptr) :: m
        type(_METHOD_), pointer :: mp

        call c_f_pointer(m, mp)
        call mp%finalize
        deallocate( mp )
    end subroutine S(c_finalize)


    function S(c_get_eigenvalues)(m, dim, which) &
        result(evp) bind(c)
        use iso_c_binding
        type(c_ptr), value :: m
        type(c_ptr) :: evp
        integer, intent(out)  :: dim(1)
        integer, value :: which
        type(_METHOD_), pointer :: mp

        call c_f_pointer(m, mp)
        select case (which)
        case (1)
           !evp = c_loc(mp%eigenvalues_r)
           evp = c_loc(mp%eigenvalues_r(mp%nfrmin))
           dim(1) = mp%nf1max-mp%nf1min+1
#ifndef _HERMITE_
        case (2)
           !evp = c_loc(mp%eigenvalues_theta)
           evp = c_loc(mp%eigenvalues_theta(mp%nfthetamin))
           dim(1) = mp%nf2max-mp%nf2min+1
#else
        case (2)
           !evp = c_loc(mp%eigenvalues_z)
           evp = c_loc(mp%eigenvalues_z(mp%nfzmin))
           dim(1) = mp%nf2max-mp%nf2min+1
        case (3)
           !evp = c_loc(mp%eigenvalues_theta)
           evp = c_loc(mp%eigenvalues_theta(mp%nfthetamin))
           dim(1) = mp%nf3max-mp%nf3min+1
#endif
        end select
    end function S(c_get_eigenvalues)

   function S(c_get_nodes)(m, dim, which) &
        result(np) bind(c)
        use iso_c_binding
        type(c_ptr), value :: m
        type(c_ptr) :: np
        integer, intent(out)  :: dim(1)
        integer, value :: which
        type(_METHOD_), pointer :: mp

        call c_f_pointer(m, mp)
        select case (which)
        case (1)
            np = c_loc(mp%g%nodes_r)
            dim(1) = mp%g%n1max-mp%g%n1min+1
        case (2)
            np = c_loc(mp%g%nodes_theta)
#ifdef _HERMITE_
            dim(1) = mp%g%n3max-mp%g%n3min+1
#else
            dim(1) = mp%g%n2max-mp%g%n2min+1
#endif
        end select
    end function S(c_get_nodes)
 

    function S(c_get_weights)(m, dim, which) &
        result(np) bind(c)
        use iso_c_binding
        type(c_ptr), value :: m
        type(c_ptr) :: np
        integer, intent(out)  :: dim(1)
        integer, value :: which
        type(_METHOD_), pointer :: mp

        call c_f_pointer(m, mp)
        select case (which)
        case (1)
            !np = c_loc(mp%g%weights_r)
            np = c_loc(mp%g%weights_r(mp%g%n1min))
            dim(1) = mp%g%n1max-mp%g%n1min+1
#ifdef _HERMITE_
        case (2)
            !np = c_loc(mp%g%weights_z)
            np = c_loc(mp%g%weights_z(mp%g%n2min))
            dim(1) = mp%g%n2max-mp%g%n2min+1
#endif
        end select
    end function S(c_get_weights)


    function S(c_get_L)(m, dim) &
        result(np) bind(c)
        use iso_c_binding
        type(c_ptr), value :: m 
        type(c_ptr)  :: np
        integer, intent(out)  :: dim(3)
        type(_METHOD_), pointer :: mp

        call c_f_pointer(m, mp)
        np = c_loc(mp%L)
        dim(1) = mp%g%nr
        dim(2) = mp%nfr
        dim(3) = mp%g%ntheta/2 + 1
    end function S(c_get_L)

#ifdef _HERMITE_
    function S(c_get_H)(m, dim) &
        result(np) bind(c)
        use iso_c_binding
        type(c_ptr), value :: m 
        type(c_ptr)  :: np
        integer, intent(out)  :: dim(2)
        type(_METHOD_), pointer :: mp

        call c_f_pointer(m, mp)
        np = c_loc(mp%H_z)
        dim(1) = mp%g%n3max-mp%g%n3min+1
        dim(2) = mp%nf3max-mp%nf3min+1
    end function S(c_get_H)
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function S(c_new_wf)(m) &
        result(this) bind(c)
        use iso_c_binding
        type(c_ptr) :: this 
        type(c_ptr), value :: m 
        type(_METHOD_), pointer :: mp
        type(_WF_), pointer :: psi

        call c_f_pointer(m, mp)
        allocate( psi )
        psi = _WF_(mp)
        this =  c_loc(psi)
    end function S(c_new_wf)

    subroutine S(c_finalize_wf)(psi)  bind(c)
        use iso_c_binding
        type(c_ptr) :: psi 
        type(_WF_), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%finalize
        deallocate( psip )
    end subroutine S(c_finalize_wf)


    function S(c_is_real_space_wf)(psi) &
        result(ans) bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        logical :: ans

        call c_f_pointer(psi, psip)
        ans = psip%is_real_space
    end function S(c_is_real_space_wf)

    subroutine S(c_to_real_space_wf)(psi) bind(c) 
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%to_real_space
    end subroutine S(c_to_real_space_wf)
  
    subroutine S(c_to_frequency_space_wf)(psi) bind(c) 
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%to_frequency_space
    end subroutine S(c_to_frequency_space_wf)

    subroutine S(c_propagate_A_wf)(psi, dt) bind(c) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(_WF_), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%propagate_A(dt)
    end subroutine S(c_propagate_A_wf)

    subroutine S(c_propagate_B_wf)(psi, dt) bind(c) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(_WF_), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%propagate_B(dt)
    end subroutine S(c_propagate_B_wf)

    subroutine S(c_propagate_C_wf)(psi, dt) bind(c) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(_WF_), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%propagate_C(dt)
    end subroutine S(c_propagate_C_wf)

    subroutine S(c_add_apply_A_wf)(this, other, coefficient) bind(c) 
        use iso_c_binding
        type(c_ptr), value :: this
        type(c_ptr), value :: other 
        _COMPLEX_OR_REAL_(kind=prec), value :: coefficient
        type(_WF_), pointer :: thisp
        type(_WF_), pointer :: otherp

        call c_f_pointer(this, thisp)
        call c_f_pointer(other, otherp)
        call thisp%add_apply_A(otherp, coefficient)
    end subroutine S(c_add_apply_A_wf)


    subroutine S(c_save_wf)(psi, filename, filename_length) bind(c) 
        use iso_c_binding
        type(c_ptr), value :: psi
        integer(c_int), intent(in), value :: filename_length
        character(kind=c_char), dimension(filename_length), intent(in):: filename
        character(len=filename_length, kind=c_char) :: fn 
        type(_WF_), pointer :: psip
        integer :: j
        call c_f_pointer(psi, psip)
        do j=1,filename_length
            fn(j:j) = filename(j)
        end do
        call psip%save(fn)
    end subroutine S(c_save_wf)

    subroutine S(c_load_wf)(psi, filename, filename_length) bind(c) 
        use iso_c_binding
        type(c_ptr), value :: psi
        integer(c_int), intent(in), value :: filename_length
        character(kind=c_char), dimension(filename_length), intent(in):: filename
        character(len=filename_length, kind=c_char) :: fn 
        type(_WF_), pointer :: psip
        integer :: j
        call c_f_pointer(psi, psip)
        do j=1,filename_length
            fn(j:j) = filename(j)
        end do
        call psip%load(fn)
    end subroutine S(c_load_wf)

    function S(c_norm_wf)(psi) &
        result(ans) bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        real(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        ans = psip%norm()
    end function S(c_norm_wf)

   function S(c_inner_product_wf)(psi, other) &
        result(ans) bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(c_ptr), value :: other 
        type(_WF_), pointer :: psip
        type(_WF_), pointer :: otherp
        _COMPLEX_OR_REAL_(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        call c_f_pointer(other, otherp)
        ans = psip%inner_product(otherp)
   end function S(c_inner_product_wf)


    function S(c_norm_in_frequency_space_wf)(psi) &
        result(ans) bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        real(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        ans = psip%norm_in_frequency_space()
    end function S(c_norm_in_frequency_space_wf)
    

   function S(c_normalize_wf)(psi) &
        result(ans) bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        real(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        call psip%normalize(ans)
    end function S(c_normalize_wf)

   subroutine S(c_scale_wf)(psi, factor) bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        _COMPLEX_OR_REAL_(kind=prec), value :: factor 

        call c_f_pointer(psi, psip)
        call psip%scale(factor)
    end subroutine S(c_scale_wf)

   subroutine S(c_axpy_wf)(this, other, factor) bind(c) 
        use iso_c_binding
        type(c_ptr), value :: this
        type(c_ptr), value :: other 
        _COMPLEX_OR_REAL_(kind=prec), value :: factor
        type(_WF_), pointer :: thisp
        type(_WF_), pointer :: otherp

        call c_f_pointer(this, thisp)
        call c_f_pointer(other, otherp)
        call thisp%axpy(otherp, factor)
    end subroutine S(c_axpy_wf)



    function S(c_get_data_wf)(psi, dim) &
        result(up) bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(c_ptr)  :: up
        integer, intent(out)  :: dim(_DIM_)
        type(_WF_), pointer :: psip

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
    end function S(c_get_data_wf)






#ifndef _REAL_ 
    subroutine S(c_rset_wf)(psi, f) &
        bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
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
    end subroutine S(c_rset_wf)
#endif    

    subroutine S(c_set_wf)(psi, f) &
        bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
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
    end subroutine S(c_set_wf)


#ifndef _REAL_ 
    subroutine S(c_rset_t_wf)(psi, f, t) &
        bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
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
    end subroutine S(c_rset_t_wf)
#endif    

    subroutine S(c_set_t_wf)(psi, f, t) &
        bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
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
    end subroutine S(c_set_t_wf)



    subroutine S(c_copy_wf)(psi, source) &
        bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(c_ptr), value :: source 
        type(_WF_), pointer :: psip, sourcep

        call c_f_pointer(psi, psip)
        call c_f_pointer(source, sourcep)
        call psip%copy(sourcep)
    end subroutine S(c_copy_wf)

    function S(c_distance_wf)(psi1, psi2) &
        result(d) bind(c)
        use iso_c_binding
        type(c_ptr), value :: psi1 
        type(c_ptr), value :: psi2 
        real(kind = prec) :: d
        type(_WF_), pointer :: psi1p, psi2p

        call c_f_pointer(psi1, psi1p)
        call c_f_pointer(psi2, psi2p)
        d =  psi1p%distance(psi2p)
    end function S(c_distance_wf)


end module _C_MODULE_ 


#endif
