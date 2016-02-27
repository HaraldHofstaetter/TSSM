#ifdef _ROTSYM_
  #define _DIM_ 1
  #ifdef _REAL_
    #define _METHOD_ bessel_rotsym_real_1d
    #define _WF_ wf_bessel_rotsym_real_1d
    #define _COMPLEX_OR_REAL_ real
    #define S(x) x ## _bessel_rotsym_real_1d
  #else
    #define _METHOD_ bessel_rotsym_1d
    #define _WF_ wf_bessel_rotsym_1d
    #define _COMPLEX_OR_REAL_ complex
    #define S(x) x ## _bessel_rotsym_1d
#endif
#else
  #define _DIM_ 2
  #ifdef _REAL_
    #define _METHOD_ fourier_bessel_real_2d
    #define _WF_ wf_fourier_bessel_real_2d
    #define S(x) x ## _fourier_bessel_real_2d
  #else
    #define _METHOD_ fourier_bessel_2d
    #define _WF_ wf_fourier_bessel_2d
    #define S(x) x ## _fourier_bessel_2d
  #endif
#endif

#ifdef _REAL_
    #define _COMPLEX_OR_REAL_ real
#else
    #define _COMPLEX_OR_REAL_ complex
#endif

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



#ifdef _QUADPRECISION_
module S(tssmq_c)
    use tssmq_base
    use tssmq_fourier_bessel
#else
module S(tssm_c)
    use tssm_base
    use tssm_fourier_bessel
#endif
    implicit none
    private

contains

#ifdef _ROTSYM_        
    function c_new(nr, r_max, boundary_conditions) &
#else
    function c_new(ntheta, nr, nfr, r_max, boundary_conditions, quadrature_formula) &
#endif    
        result(this) bind(c, name=SC(new))
        use iso_c_binding
        type(c_ptr) :: this 
        integer, value ::nr 
        real(kind=prec),  value :: r_max 
        integer, value :: boundary_conditions
#ifndef _ROTSYM_        
        integer, value :: ntheta
        integer, value ::nfr 
        integer, value :: quadrature_formula
#endif        

        type(_METHOD_), pointer :: meth
        allocate( meth )
#ifdef _ROTSYM_        
        meth = _METHOD_( nr, r_max, boundary_conditions)
#else
        meth = _METHOD_( ntheta, nr, nfr, r_max, boundary_conditions, quadrature_formula)
#endif        
        this =  c_loc(meth)
    end function c_new


#ifndef _ROTSYM_        
    function c_new_from_file(filename, filename_length) &
        result(this) bind(c, name=SC(new_from_file))
        use iso_c_binding
        type(c_ptr) :: this 
        integer(c_int), intent(in), value :: filename_length
        character(kind=c_char), dimension(filename_length), intent(in):: filename
        character(len=filename_length, kind=c_char) :: fn 
        integer :: j
        type(_METHOD_), pointer :: meth

        allocate( meth )
        do j=1,filename_length
            fn(j:j) = filename(j)
        end do
        meth = _METHOD_(fn)
        this =  c_loc(meth)
    end function c_new_from_file
#endif    

    subroutine c_finalize(m)  bind(c, name=SC(finalize))
        use iso_c_binding
        type(c_ptr) :: m
        type(_METHOD_), pointer :: mp

        call c_f_pointer(m, mp)
        call mp%finalize
        deallocate( mp )
    end subroutine c_finalize

   subroutine c_set_propagate_time_together_with_A(m, flag) &
        bind(c, name=SC(set_propagate_time_together_with_A)) 
        use iso_c_binding
        type(c_ptr), value :: m 
        logical, value :: flag 
        type(_METHOD_), pointer :: mp

        call c_f_pointer(m, mp)
        mp%propagate_time_together_with_A = flag 
    end subroutine c_set_propagate_time_together_with_A

    function c_get_propagate_time_together_with_A(m) &
        result(flag) bind(c, name=SC(get_propagate_time_together_with_A))
        use iso_c_binding
        type(c_ptr), value :: m 
        type(_METHOD_), pointer :: mp
        logical :: flag 

        call c_f_pointer(m, mp)
        flag = mp%propagate_time_together_with_A
    end function c_get_propagate_time_together_with_A

#ifndef _ROTSYM_        
    subroutine c_save_method(m, filename, filename_length) &
        bind(c, name=SC(save)) 
        use iso_c_binding
        type(c_ptr), value :: m
        integer(c_int), intent(in), value :: filename_length
        character(kind=c_char), dimension(filename_length), intent(in):: filename
        character(len=filename_length, kind=c_char) :: fn 
        type(_METHOD_), pointer :: mp
        integer :: j
        call c_f_pointer(m, mp)
        do j=1,filename_length
            fn(j:j) = filename(j)
        end do
        call mp%save(fn)
    end subroutine c_save_method
#endif    
    

   function c_get_nr(m) &
        result(nr) bind(c, name=SC(get_nr))
        use iso_c_binding
        type(c_ptr), value :: m
        integer :: nr
        type(_METHOD_), pointer :: mp
        call c_f_pointer(m, mp)
#ifdef _ROTSYM_
        nr = mp%g%nx
#else
        nr = mp%g%nr
#endif        
    end function c_get_nr       

    function c_get_rmax(m) &
        result(rmax) bind(c, name=SC(get_rmax))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: rmax
        type(_METHOD_), pointer :: mp
        call c_f_pointer(m, mp)
        rmax = mp%rmax
    end function c_get_rmax   
    
#ifndef _ROTSYM_
    function c_get_ntheta(m) &
        result(ntheta) bind(c, name=SC(get_ntheta))
        use iso_c_binding
        type(c_ptr), value :: m
        integer :: ntheta
        type(_METHOD_), pointer :: mp
        call c_f_pointer(m, mp)
        ntheta = mp%g%ntheta
    end function c_get_ntheta   
    
    function c_get_nfr(m) &
        result(nfr) bind(c, name=SC(get_nfr))
        use iso_c_binding
        type(c_ptr), value :: m
        integer :: nfr
        type(_METHOD_), pointer :: mp
        call c_f_pointer(m, mp)
        nfr = mp%nfr
    end function c_get_nfr   
#endif    

    function c_get_eigenvalues(m, dim) &
        result(evp) bind(c, name=SC(get_eigenvalues))
        use iso_c_binding
        type(c_ptr), value :: m
        type(c_ptr) :: evp
        integer, intent(out)  :: dim(_DIM_)
        type(_METHOD_), pointer :: mp

        call c_f_pointer(m, mp)
#ifdef _ROTSYM_
        evp = c_loc(mp%eigenvalues1(mp%nf1min))
#else
        !evp = c_loc(mp%eigenvalues_r_theta)
        evp = c_loc(mp%eigenvalues_r_theta(mp%nf1min,mp%nf2min))
#endif        
        dim(1) = mp%nf1max-mp%nf1min+1
#ifndef _ROTSYM_
        dim(2) = mp%nf2max-mp%nf2min+1
#endif        
    end function c_get_eigenvalues

    function c_get_nodes(m, dim, which) &
        result(np) bind(c, name=SC(get_nodes))
        use iso_c_binding
        type(c_ptr), value :: m
        type(c_ptr) :: np
        integer, intent(out)  :: dim(_DIM_)
        integer, value :: which
        type(_METHOD_), pointer :: mp

        call c_f_pointer(m, mp)
        select case (which)
        case (1)
#ifdef _ROTSYM_
            np = c_loc(mp%g%nodes_x)
#else
            np = c_loc(mp%g%nodes_r)
#endif            
            dim(1) = mp%g%n1max-mp%g%n1min+1
#ifndef _ROTSYM_
        case (2)
            np = c_loc(mp%g%nodes_theta)
            dim(1) = mp%g%n2max-mp%g%n2min+1
#endif            
        end select
    end function c_get_nodes
 

    function c_get_weights(m, dim) &
        result(np) bind(c, name=SC(get_weights))
        use iso_c_binding
        type(c_ptr), value :: m
        type(c_ptr) :: np
        integer, intent(out)  :: dim(1)
        type(_METHOD_), pointer :: mp

        call c_f_pointer(m, mp)
#ifdef _ROTSYM_
        np = c_loc(mp%g%weights_x(mp%g%n1min))
#else
        !np = c_loc(mp%g%weights_r)
        np = c_loc(mp%g%weights_r(mp%g%n1min))
#endif        
        dim(1) = mp%g%n1max-mp%g%n1min+1
    end function c_get_weights


    function c_get_L(m, dim) &
        result(np) bind(c, name=SC(get_L))
        use iso_c_binding
        type(c_ptr), value :: m 
        type(c_ptr)  :: np
        integer, intent(out)  :: dim(_DIM_+1)
        type(_METHOD_), pointer :: mp

        call c_f_pointer(m, mp)
#ifdef _ROTSYM_
        np = c_loc(mp%H1)
        dim(1) = mp%g%nx
        dim(2) = mp%nf1
#else
        np = c_loc(mp%L)
        dim(1) = mp%g%nr
        dim(2) = mp%nfr
        dim(3) = mp%g%ntheta/2 + 1
#endif        
    end function c_get_L

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function c_new_wf(m) &
        result(this) bind(c, name=SC(new_wf))
        use iso_c_binding
        type(c_ptr) :: this 
        type(c_ptr), value :: m 
        type(_METHOD_), pointer :: mp
        type(_WF_), pointer :: psi

        call c_f_pointer(m, mp)
        allocate( psi )
        psi = _WF_(mp)
        this =  c_loc(psi)
    end function c_new_wf

    subroutine c_finalize_wf(psi) &
        bind(c, name=SC(finalize_wf))
        use iso_c_binding
        type(c_ptr) :: psi 
        type(_WF_), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%finalize
        deallocate( psip )
    end subroutine c_finalize_wf


    function c_is_real_space_wf(psi) &
        result(ans) bind(c, name=SC(is_real_space_wf))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        logical :: ans

        call c_f_pointer(psi, psip)
        ans = psip%is_real_space
    end function c_is_real_space_wf

    subroutine c_to_real_space_wf(psi) &
        bind(c, name=SC(to_real_space_wf)) 
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%to_real_space
    end subroutine c_to_real_space_wf
  
    subroutine c_to_frequency_space_wf(psi) &
        bind(c, name=SC(to_frequency_space_wf)) 
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%to_frequency_space
    end subroutine c_to_frequency_space_wf

    subroutine c_set_time_wf(psi, t) &
        bind(c, name=SC(set_time_wf)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: t
        type(_WF_), pointer :: psip

        call c_f_pointer(psi, psip)
        psip%time = t
    end subroutine c_set_time_wf

    function c_get_time_wf(psi) &
        result(ans) bind(c, name=SC(get_time_wf))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        _COMPLEX_OR_REAL_(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        ans = psip%time
    end function c_get_time_wf

    subroutine c_propagate_time_wf(psi, dt) &
        bind(c, name=SC(propagate_time_wf)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(_WF_), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%propagate_time(dt)
    end subroutine c_propagate_time_wf
    

    subroutine c_propagate_A_wf(psi, dt) &
        bind(c, name=SC(propagate_A_wf)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(_WF_), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%propagate_A(dt)
    end subroutine c_propagate_A_wf

    subroutine c_propagate_A_derivative_wf(psi, wf, dt) & 
        bind(c, name=SC(propagate_A_derivative_wf)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        type(c_ptr), value :: wf
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(_WF_), pointer :: psip
        type(_WF_), pointer :: wfp 

        call c_f_pointer(psi, psip)
        call c_f_pointer(wf, wfp)
        call psip%propagate_A_derivative(wfp, dt)
    end subroutine c_propagate_A_derivative_wf
    

    subroutine c_propagate_B_wf(psi, dt) &
        bind(c, name=SC(propagate_B_wf)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(_WF_), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%propagate_B(dt)
    end subroutine c_propagate_B_wf

    subroutine c_propagate_C_wf(psi, dt) &
        bind(c, name=SC(propagate_C_wf)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(_WF_), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%propagate_C(dt)
    end subroutine c_propagate_C_wf

    subroutine c_add_apply_A_wf(this, other, coefficient) &
        bind(c, name=SC(add_apply_A_wf)) 
        use iso_c_binding
        type(c_ptr), value :: this
        type(c_ptr), value :: other 
        _COMPLEX_OR_REAL_(kind=prec), value :: coefficient
        type(_WF_), pointer :: thisp
        type(_WF_), pointer :: otherp

        call c_f_pointer(this, thisp)
        call c_f_pointer(other, otherp)
        call thisp%add_apply_A(otherp, coefficient)
    end subroutine c_add_apply_A_wf


    subroutine c_save_wf(psi, filename, filename_length) &
        bind(c, name=SC(save_wf)) 
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
    end subroutine c_save_wf

    subroutine c_load_wf(psi, filename, filename_length) &
        bind(c, name=SC(load_wf)) 
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
    end subroutine c_load_wf

#ifndef _ROTSYM_
   function c_evaluate_wf(psi, x, y) &
        result(ans) bind(c, name=SC(evaluate_wf))
        use iso_c_binding
        type(c_ptr), value :: psi 
        real(kind=prec), value :: x 
        real(kind=prec), value :: y 
        type(_WF_), pointer :: psip
        _COMPLEX_OR_REAL_(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        ans = psip%evaluate(x, y)
    end function c_evaluate_wf
#endif    
    

    function c_norm_wf(psi) &
        result(ans) bind(c, name=SC(norm_wf))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        real(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        ans = psip%norm()
    end function c_norm_wf

    function c_norm_in_frequency_space_wf(psi) &
        result(ans) bind(c, name=SC(norm_in_frequency_space_wf))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        real(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        ans = psip%norm_in_frequency_space()
    end function c_norm_in_frequency_space_wf

#ifndef _ROTSYM_
    function c_inner_product_wf(psi, other) &
        result(ans) bind(c, name=SC(inner_product_wf))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(c_ptr), value :: other 
        type(_WF_), pointer :: psip
        type(_WF_), pointer :: otherp
        _COMPLEX_OR_REAL_(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        call c_f_pointer(other, otherp)
        ans = psip%inner_product(otherp)
   end function c_inner_product_wf
#endif    

   function c_normalize_wf(psi) &
        result(ans) bind(c, name=SC(normalize_wf))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        real(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        call psip%normalize(ans)
    end function c_normalize_wf

   subroutine c_scale_wf(psi, factor) bind(c, name=SC(scale_wf))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        _COMPLEX_OR_REAL_(kind=prec), value :: factor 

        call c_f_pointer(psi, psip)
        call psip%scale(factor)
    end subroutine c_scale_wf

   subroutine c_axpy_wf(this, other, factor) bind(c, name=SC(axpy_wf)) 
        use iso_c_binding
        type(c_ptr), value :: this
        type(c_ptr), value :: other 
        _COMPLEX_OR_REAL_(kind=prec), value :: factor
        type(_WF_), pointer :: thisp
        type(_WF_), pointer :: otherp

        call c_f_pointer(this, thisp)
        call c_f_pointer(other, otherp)
        call thisp%axpy(otherp, factor)
    end subroutine c_axpy_wf



    function c_get_data_wf(psi, dim) &
        result(up) bind(c, name=SC(get_data_wf))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(c_ptr)  :: up
        integer, intent(out)  :: dim(_DIM_)
        type(_WF_), pointer :: psip

        call c_f_pointer(psi, psip)
        !up = c_loc(psip%u)
        !up = c_loc(psip%u(psip%m%g%m1min, psip%m%g%m2min))
        up = c_loc(psip%up(1))
        dim(1) = psip%m%g%m1max-psip%m%g%m1min+1
#ifndef _ROTSYM_
        dim(2) = psip%m%g%m2max-psip%m%g%m2min+1
#endif        
    end function c_get_data_wf






#ifndef _REAL_ 
#ifdef _ROTSYM_
    subroutine c_rset_wf(psi, f) &
        bind(c, name=SC(rset_wf))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        interface 
           function f(r) bind(c)
               import prec
               real(kind=prec), value :: r
               real(kind=prec) :: f 
           end function f
        end interface 

        call c_f_pointer(psi, psip)
        call psip%rset(ff)
    contains
        function ff(r)
            real(kind=prec), intent(in) :: r
            real(kind=prec) :: ff 
            ff = f(r)
        end function ff
    end subroutine c_rset_wf
#else
    subroutine c_rset_wf(psi, f) &
        bind(c, name=SC(rset_wf))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        interface 
           function f(x, y) bind(c)
               import prec
               real(kind=prec), value :: x
               real(kind=prec), value :: y
               real(kind=prec) :: f 
           end function f
        end interface 

        call c_f_pointer(psi, psip)
        call psip%rset(ff)
    contains
        function ff(x, y)
            real(kind=prec), intent(in) :: x
            real(kind=prec), intent(in) :: y
            real(kind=prec) :: ff 
            ff = f(x, y)
        end function ff
    end subroutine c_rset_wf
#endif    
#endif    

#ifdef _ROTSYM_
    subroutine c_set_wf(psi, f) &
        bind(c, name=SC(set_wf))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        interface 
           function f(r) bind(c)
               import prec
               real(kind=prec), value :: r
               _COMPLEX_OR_REAL_(kind=prec) :: f 
           end function f
        end interface 

        call c_f_pointer(psi, psip)
        call psip%set(ff)
    contains
        function ff(r)
            real(kind=prec), intent(in) :: r
            _COMPLEX_OR_REAL_(kind=prec) :: ff 
            ff = f(r)
        end function ff
    end subroutine c_set_wf

#else
    subroutine c_set_wf(psi, f) &
        bind(c, name=SC(set_wf))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        interface 
           function f(x, y) bind(c)
               import prec
               real(kind=prec), value :: x
               real(kind=prec), value :: y
               _COMPLEX_OR_REAL_(kind=prec) :: f 
           end function f
        end interface 

        call c_f_pointer(psi, psip)
        call psip%set(ff)
    contains
        function ff(x, y)
            real(kind=prec), intent(in) :: x
            real(kind=prec), intent(in) :: y
            _COMPLEX_OR_REAL_(kind=prec) :: ff 
            ff = f(x, y)
        end function ff
    end subroutine c_set_wf
#endif    


#ifndef _REAL_ 
#ifdef _ROTSYM_
    subroutine c_rset_t_wf(psi, f, t) &
        bind(c, name=SC(rset_t_wf))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        real(kind=prec), value :: t
        interface 
           function f(r, t) bind(c)
               import prec
               real(kind=prec), value :: r
               real(kind=prec) :: f 
               real(kind=prec), value :: t
           end function f
        end interface 

        call c_f_pointer(psi, psip)
        call psip%rset_t(ff, t)
    contains
        function ff(r, t)
            real(kind=prec), intent(in) :: r
            real(kind=prec), intent(in) :: t
            real(kind=prec) :: ff 
            ff = f(r, t)
        end function ff
    end subroutine c_rset_t_wf
#else
    subroutine c_rset_t_wf(psi, f, t) &
        bind(c, name=SC(rset_t_wf))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        real(kind=prec), value :: t
        interface 
           function f(x, y, t) bind(c)
               import prec
               real(kind=prec), value :: x
               real(kind=prec), value :: y
               real(kind=prec) :: f 
               real(kind=prec), value :: t
           end function f
        end interface 

        call c_f_pointer(psi, psip)
        call psip%rset_t(ff, t)
    contains
        function ff(x, y, t)
            real(kind=prec), intent(in) :: x
            real(kind=prec), intent(in) :: y
            real(kind=prec), intent(in) :: t
            real(kind=prec) :: ff 
            ff = f(x, y, t)
        end function ff
    end subroutine c_rset_t_wf
#endif    
#endif    

#ifdef _ROTSYM_
    subroutine c_set_t_wf(psi, f, t) &
        bind(c, name=SC(set_t_wf))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        real(kind=prec), value :: t
        interface 
           function f(r, t) bind(c)
               import prec
               real(kind=prec), value :: r
               _COMPLEX_OR_REAL_(kind=prec) :: f 
               real(kind=prec), value :: t
           end function f
        end interface 

        call c_f_pointer(psi, psip)
        call psip%set_t(ff, t)
    contains
        function ff(r, t)
            real(kind=prec), intent(in) :: r
            real(kind=prec), intent(in) :: t
            _COMPLEX_OR_REAL_(kind=prec) :: ff 
            ff = f(r, t)
        end function ff
    end subroutine c_set_t_wf
#else
    subroutine c_set_t_wf(psi, f, t) &
        bind(c, name=SC(set_t_wf))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(_WF_), pointer :: psip
        real(kind=prec), value :: t
        interface 
           function f(x, y, t) bind(c)
               import prec
               real(kind=prec), value :: x
               real(kind=prec), value :: y
               _COMPLEX_OR_REAL_(kind=prec) :: f 
               real(kind=prec), value :: t
           end function f
        end interface 

        call c_f_pointer(psi, psip)
        call psip%set_t(ff, t)
    contains
        function ff(x, y, t)
            real(kind=prec), intent(in) :: x
            real(kind=prec), intent(in) :: y
            real(kind=prec), intent(in) :: t
            _COMPLEX_OR_REAL_(kind=prec) :: ff 
            ff = f(x, y, t)
        end function ff
    end subroutine c_set_t_wf
#endif



    subroutine c_copy_wf(psi, source) &
        bind(c, name=SC(copy_wf))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(c_ptr), value :: source 
        type(_WF_), pointer :: psip, sourcep

        call c_f_pointer(psi, psip)
        call c_f_pointer(source, sourcep)
        call psip%copy(sourcep)
    end subroutine c_copy_wf

    function c_distance_wf(psi1, psi2) &
        result(d) bind(c, name=SC(distance_wf))
        use iso_c_binding
        type(c_ptr), value :: psi1 
        type(c_ptr), value :: psi2 
        real(kind = prec) :: d
        type(_WF_), pointer :: psi1p, psi2p

        call c_f_pointer(psi1, psi1p)
        call c_f_pointer(psi2, psi2p)
        d =  psi1p%distance(psi2p)
    end function c_distance_wf


#ifdef _QUADPRECISION_
end module S(tssmq_c)
#else
end module S(tssm_c)
#endif
