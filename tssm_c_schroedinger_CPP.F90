#ifdef _HERMITE_
#ifdef _REAL_
#define S0(x,y)  x ## _hermite_real_ ## y ## d 
#else
#define S0(x,y)  x ## _hermite_ ## y ## d 
#endif
#else
#ifdef _REAL_
#define S0(x,y)  x ## _real_ ## y ## d 
#else
#define S0(x,y)  x ## _ ## y ## d 
#endif
#endif
#define S1(x,y) S0(x,y)
#define S(x) S1(x,_DIM_)

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
#define SB0(x,y)  x ## _real_ ## y ## d 
#else
#define SB0(x,y)  x ## _ ## y ## d 
#endif
#define SB1(x,y) SB0(x,y)
#define SB(x) SB1(x,_DIM_)

#ifdef _REAL_
#define _WAVE_FUNCTION_ real_wave_function
#define _COMPLEX_OR_REAL_ real
#else
#define _WAVE_FUNCTION_ wave_function
#define _COMPLEX_OR_REAL_ complex
#endif



#ifdef _QUADPRECISION_
module S(tssmq_c_schroedinger)
    use tssmq_schroedinger
#else
module S(tssm_c_schroedinger)
    use tssm_schroedinger
#endif    
    implicit none

contains


#ifdef _HERMITE_

   function S(c_new_schroedinger)(nx, omega_x, &
#if(_DIM_>=2)
                                  ny, omega_y, &
#endif
#if(_DIM_>=3)
                                  nz, omega_z, &
#endif
                              hbar, mass, potential, with_potential, cubic_coupling) &
        result(this) bind(c, name=SC(new_schroedinger))
        use iso_c_binding
        type(c_ptr) :: this 
        real(kind=prec),  value :: omega_x 
        integer, value :: nx
#if(_DIM_>=2)
        real(kind=prec),  value :: omega_y 
        integer, value :: ny
#endif
#if(_DIM_>=3)
        real(kind=prec),  value :: omega_z 
        integer, value :: nz
#endif
        real(kind=prec), value  :: hbar 
        real(kind=prec), value  :: mass 
        real(kind=prec), external, bind(c) :: potential 
        logical, value  :: with_potential
        real(kind=prec), value  :: cubic_coupling 

        type(S(schroedinger)), pointer :: m
        allocate( m )
        m = S(schroedinger)(nx, omega_x, &
#if(_DIM_>=2)
                            ny, omega_y, &
#endif
#if(_DIM_>=3)
                            nz, omega_z, &
#endif
                       hbar, mass, cubic_coupling=cubic_coupling) 
        this =  c_loc(m)
        if (with_potential) then
           call S(c_set_potential_schroedinger)(this, potential)
        end if
    end function S(c_new_schroedinger)

#else

    function S(c_new_schroedinger)(nx, xmin, xmax, &
#if(_DIM_>=2)
                              ny, ymin, ymax, &
#endif
#if(_DIM_>=3)
                              nz, zmin, zmax, &
#endif
                              hbar, mass, potential, with_potential, cubic_coupling, boundary_conditions) &
        result(this) bind(c, name=SC(new_schroedinger))
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
        real(kind=prec), value  :: hbar 
        real(kind=prec), value  :: mass 
        real(kind=prec), external, bind(c) :: potential 
        logical, value  :: with_potential
        real(kind=prec), value  :: cubic_coupling 
        integer, value :: boundary_conditions

        type(S(schroedinger)), pointer :: m
        allocate( m )
        m = S(schroedinger)(nx, xmin, xmax, &
#if(_DIM_>=2)
                       ny, ymin, ymax, &
#endif
#if(_DIM_>=3)
                       nz, zmin, zmax, &
#endif
                       hbar, mass, cubic_coupling=cubic_coupling, boundary_conditions=boundary_conditions) 
        this =  c_loc(m)
        if (with_potential) then
            call S(c_set_potential_schroedinger)(this, potential)          
        end if
    end function S(c_new_schroedinger)
#endif


    subroutine S(c_finalize_schroedinger)(m) &
        bind(c, name=SC(finalize_schroedinger))
        use iso_c_binding
        type(c_ptr) :: m
        type(S(schroedinger)), pointer :: mp

        call c_f_pointer(m, mp)
        call mp%finalize
        deallocate( mp )
    end subroutine S(c_finalize_schroedinger)



    function S(c_new_wf_schroedinger)(m) &
        result(this) bind(c, name=SC(new_wf_schroedinger))
        use iso_c_binding
        type(c_ptr) :: this 
        type(c_ptr), value :: m 
        type(S(schroedinger)), pointer :: mp
        type(S(wf_schroedinger)), pointer :: psi

        call c_f_pointer(m, mp)
        allocate( psi )
        psi = S(wf_schroedinger)(mp)
        this =  c_loc(psi)
    end function S(c_new_wf_schroedinger)

    subroutine S(c_finalize_wfs)(psi)  &
        bind(c, name=SC(finalize_wf_schroedinger))
        use iso_c_binding
        type(c_ptr) :: psi 
        type(S(wf_schroedinger)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%finalize
        deallocate( psip )
    end subroutine S(c_finalize_wfs)


    function S(c_is_real_space_wfs)(psi) &
        result(ans) bind(c, name=SC(is_real_space_wf_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_schroedinger)), pointer :: psip
        logical :: ans

        call c_f_pointer(psi, psip)
        ans = psip%is_real_space
    end function S(c_is_real_space_wfs)

    subroutine S(c_to_real_space_wfs)(psi) &
        bind(c, name=SC(to_real_space_wf_schroedinger)) 
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_schroedinger)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%to_real_space
    end subroutine S(c_to_real_space_wfs)
  
    subroutine S(c_to_frequency_space_wfs)(psi) & 
        bind(c, name=SC(to_frequency_space_wf_schroedinger)) 
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_schroedinger)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%to_frequency_space
    end subroutine S(c_to_frequency_space_wfs)

    subroutine S(c_propagate_A_wfs)(psi, dt) & 
         bind(c, name=SC(propagate_A_wf_schroedinger)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(S(wf_schroedinger)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%propagate_A(dt)
    end subroutine S(c_propagate_A_wfs)

    subroutine S(c_propagate_B_wfs)(psi, dt) & 
        bind(c, name=SC(propagate_B_wf_schroedinger)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(S(wf_schroedinger)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%propagate_B(dt)
    end subroutine S(c_propagate_B_wfs)

    subroutine S(c_propagate_C_wfs)(psi, dt) & 
        bind(c, name=SC(propagate_C_wf_schroedinger)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(S(wf_schroedinger)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%propagate_C(dt)
    end subroutine S(c_propagate_C_wfs)
    
    subroutine S(c_add_apply_A_wfs)(this, other, coefficient) &
        bind(c, name=SC(add_apply_A_wf_schroedinger)) 
        use iso_c_binding
        type(c_ptr), value :: this
        type(c_ptr), value :: other 
        _COMPLEX_OR_REAL_(kind=prec), value :: coefficient
        type(S(wf_schroedinger)), pointer :: thisp
        type(S(wf_schroedinger)), pointer :: otherp

        call c_f_pointer(this, thisp)
        call c_f_pointer(other, otherp)
        call thisp%add_apply_A(otherp, coefficient)
    end subroutine S(c_add_apply_A_wfs)


    subroutine S(c_save_wfs)(psi, filename, filename_length) &
        bind(c, name=SC(save_wf_schroedinger)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        integer(c_int), intent(in), value :: filename_length
        character(kind=c_char), dimension(filename_length), intent(in):: filename
        character(len=filename_length, kind=c_char) :: fn 
        type(S(wf_schroedinger)), pointer :: psip
        integer :: j
        call c_f_pointer(psi, psip)
        do j=1,filename_length
            fn(j:j) = filename(j)
        end do
        call psip%save(fn)
    end subroutine S(c_save_wfs)

    subroutine S(c_load_wfs)(psi, filename, filename_length) &
        bind(c, name=SC(load_wf_schroedinger)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        integer(c_int), intent(in), value :: filename_length
        character(kind=c_char), dimension(filename_length), intent(in):: filename
        character(len=filename_length, kind=c_char) :: fn 
        type(S(wf_schroedinger)), pointer :: psip
        integer :: j
        call c_f_pointer(psi, psip)
        do j=1,filename_length
            fn(j:j) = filename(j)
        end do
        call psip%load(fn)
    end subroutine S(c_load_wfs)

    function S(c_norm_wfs)(psi) &
        result(ans) bind(c, name=SC(norm_wf_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_schroedinger)), pointer :: psip
        real(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        ans = psip%norm()
    end function S(c_norm_wfs)

    function S(c_norm_in_frequency_space_wfs)(psi) &
        result(ans) bind(c, name=SC(norm_in_frequency_space_wf_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_schroedinger)), pointer :: psip
        real(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        ans = psip%norm_in_frequency_space()
    end function S(c_norm_in_frequency_space_wfs)
    

   function S(c_normalize_wfs)(psi) &
        result(ans) bind(c, name=SC(normalize_wf_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_schroedinger)), pointer :: psip
        real(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        call psip%normalize(ans)
    end function S(c_normalize_wfs)

   subroutine S(c_scale_wfs)(psi, factor) & 
        bind(c, name=SC(scale_wf_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_schroedinger)), pointer :: psip
        _COMPLEX_OR_REAL_(kind=prec), value :: factor 

        call c_f_pointer(psi, psip)
        call psip%scale(factor)
    end subroutine S(c_scale_wfs)

   subroutine S(c_axpy_wfs)(this, other, factor) &
        bind(c, name=SC(axpy_wf_schroedinger)) 
        use iso_c_binding
        type(c_ptr), value :: this
        type(c_ptr), value :: other 
        _COMPLEX_OR_REAL_(kind=prec), value :: factor
        type(S(wf_schroedinger)), pointer :: thisp
        type(S(wf_schroedinger)), pointer :: otherp

        call c_f_pointer(this, thisp)
        call c_f_pointer(other, otherp)
        call thisp%axpy(otherp, factor)
    end subroutine S(c_axpy_wfs)

    function S(c_get_data_wfs)(psi, dim) &
        result(up) bind(c, name=SC(get_data_wf_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(c_ptr)  :: up
        integer, intent(out)  :: dim(_DIM_)
        type(S(wf_schroedinger)), pointer :: psip

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
    end function S(c_get_data_wfs)


    function S(c_get_eigenvalues_schroedinger)(m, dim, which) &
        result(evp) bind(c, name=SC(get_eigenvalues_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m
        type(c_ptr) :: evp
        integer, intent(out)  :: dim(1)
        integer, value :: which
        type(S(schroedinger)), pointer :: mp

        call c_f_pointer(m, mp)
        select case (which)
        case (1)
           !evp = c_loc(mp%eigenvalues1)
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
    end function S(c_get_eigenvalues_schroedinger)


    function S(c_get_nodes_schroedinger)(m, dim, which) &
        result(np) bind(c, name=SC(get_nodes_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m
        type(c_ptr) :: np
        integer, intent(out)  :: dim(1)
        integer, value :: which
        type(S(schroedinger)), pointer :: mp

        call c_f_pointer(m, mp)
        select case (which)
        case (1)
            np = c_loc(mp%g%nodes_x)
            dim(1) = mp%g%n1max-mp%g%n1min+1
#if(_DIM_>=2)        
        case (2)
            np = c_loc(mp%g%nodes_y)
            dim(1) = mp%g%n2max-mp%g%n2min+1
#endif
#if(_DIM_>=3)        
        case (3)
            np = c_loc(mp%g%nodes_z)
            dim(1) = mp%g%n3max-mp%g%n3min+1
#endif
        end select
    end function S(c_get_nodes_schroedinger)


!XXXXXXXXXXXXXXXXXXXXXXXXXXXX
   function S(c_get_hbar_schroedinger)(m) &
        result(hbar) bind(c, name=SC(get_hbar_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: hbar
        type(S(schroedinger)), pointer :: mp
        call c_f_pointer(m, mp)
        hbar = mp%hbar
    end function S(c_get_hbar_schroedinger)   

    function S(c_get_mass_schroedinger)(m) &
        result(mass) bind(c, name=SC(get_mass_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: mass
        type(S(schroedinger)), pointer :: mp
        call c_f_pointer(m, mp)
        mass = mp%mass
    end function S(c_get_mass_schroedinger)   

    function S(c_get_cubic_coupling_schroedinger)(m) &
        result(cubic_coupling) bind(c, name=SC(get_cubic_coupling_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: cubic_coupling
        type(S(schroedinger)), pointer :: mp
        call c_f_pointer(m, mp)
        cubic_coupling = mp%cubic_coupling
    end function S(c_get_cubic_coupling_schroedinger)   

    function S(c_get_nx_schroedinger)(m) &
        result(nx) bind(c, name=SC(get_nx_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m
        integer :: nx
        type(S(schroedinger)), pointer :: mp
        call c_f_pointer(m, mp)
        nx = mp%g%nx
    end function S(c_get_nx_schroedinger)   

#if(_DIM_>=2)
    function S(c_get_ny_schroedinger)(m) &
        result(ny) bind(c, name=SC(get_ny_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m
        integer :: ny
        type(S(schroedinger)), pointer :: mp
        call c_f_pointer(m, mp)
        ny = mp%g%ny
    end function S(c_get_ny_schroedinger)   
#endif

#if(_DIM_>=3)
    function S(c_get_nz_schroedinger)(m) &
        result(nz) bind(c, name=SC(get_nz_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m
        integer :: nz
        type(S(schroedinger)), pointer :: mp
        call c_f_pointer(m, mp)
        nz = mp%g%nz
    end function S(c_get_nz_schroedinger)   
#endif


#ifdef _HERMITE_
   function S(c_get_omega_x_schroedinger)(m) &
        result(omega_x) bind(c, name=SC(get_omega_x_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: omega_x
        type(S(schroedinger)), pointer :: mp
        call c_f_pointer(m, mp)
        omega_x = mp%omega_x
    end function S(c_get_omega_x_schroedinger)   

#if(_DIM_>=2)
   function S(c_get_omega_y_schroedinger)(m) &
        result(omega_y) bind(c, name=SC(get_omega_y_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: omega_y
        type(S(schroedinger)), pointer :: mp
        call c_f_pointer(m, mp)
        omega_y = mp%omega_y
    end function S(c_get_omega_y_schroedinger)   
#endif

#if(_DIM_>=3)
   function S(c_get_omega_z_schroedinger)(m) &
        result(omega_z) bind(c, name=SC(get_omega_z_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: omega_z
        type(S(schroedinger)), pointer :: mp
        call c_f_pointer(m, mp)
        omega_z = mp%omega_z
    end function S(c_get_omega_z_schroedinger)   
#endif

#else
   function S(c_get_xmin_schroedinger)(m) &
        result(xmin) bind(c, name=SC(get_xmin_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: xmin
        type(S(schroedinger)), pointer :: mp
        call c_f_pointer(m, mp)
        xmin = mp%g%xmin
    end function S(c_get_xmin_schroedinger)   

    function S(c_get_xmax_schroedinger)(m) &
        result(xmax) bind(c, name=SC(get_xmax_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: xmax
        type(S(schroedinger)), pointer :: mp
        call c_f_pointer(m, mp)
        xmax = mp%g%xmax
    end function S(c_get_xmax_schroedinger)     

#if(_DIM_>=2)
   function S(c_get_ymin_schroedinger)(m) &
        result(ymin) bind(c, name=SC(get_ymin_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: ymin
        type(S(schroedinger)), pointer :: mp
        call c_f_pointer(m, mp)
        ymin = mp%g%ymin
    end function S(c_get_ymin_schroedinger)   

    function S(c_get_ymax_schroedinger)(m) &
        result(ymax) bind(c, name=SC(get_ymax_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: ymax
        type(S(schroedinger)), pointer :: mp
        call c_f_pointer(m, mp)
        ymax = mp%g%ymax
    end function S(c_get_ymax_schroedinger)     
#endif

#if(_DIM_>=3)
   function S(c_get_zmin_schroedinger)(m) &
        result(zmin) bind(c, name=SC(get_zmin_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: zmin
        type(S(schroedinger)), pointer :: mp
        call c_f_pointer(m, mp)
        zmin = mp%g%zmin
    end function S(c_get_zmin_schroedinger)   

    function S(c_get_zmax_schroedinger)(m) &
        result(zmax) bind(c, name=SC(get_zmax_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m
        real(kind=prec) :: zmax
        type(S(schroedinger)), pointer :: mp
        call c_f_pointer(m, mp)
        zmax = mp%g%zmax
    end function S(c_get_zmax_schroedinger)     
#endif

#endif
!XXXXXXXXXXXXXXXXXXXXXXXXXXXX

#ifdef _HERMITE_

    function S(c_get_weights_schroedinger)(m, dim, which) &
        result(np) bind(c, name=SC(get_weights_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m
        type(c_ptr) :: np
        integer, intent(out)  :: dim(1)
        integer, value :: which
        type(S(schroedinger)), pointer :: mp

        call c_f_pointer(m, mp)
        select case (which)
        case (1)
            !np = c_loc(mp%g%weights_x)
            np = c_loc(mp%g%weights_x(mp%g%n1min))
            dim(1) = mp%g%n1max-mp%g%n1min+1
#if(_DIM_>=2)        
        case (2)
            !np = c_loc(mp%g%weights_y)
            np = c_loc(mp%g%weights_y(mp%g%n2min))
            dim(1) = mp%g%n2max-mp%g%n2min+1
#endif
#if(_DIM_>=3)        
        case (3)
            !np = c_loc(mp%g%weights_z)
            np = c_loc(mp%g%weights_z(mp%g%n3min))
            dim(1) = mp%g%n3max-mp%g%n3min+1
#endif
        end select
    end function S(c_get_weights_schroedinger)


    function S(c_get_H_schroedinger)(m, dim, which) &
        result(np) bind(c, name=SC(get_H_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m 
        type(c_ptr)  :: np
        integer, intent(out)  :: dim(2)
        integer, value :: which
        type(S(schroedinger)), pointer :: mp

        call c_f_pointer(m, mp)
        select case (which)
        case (1)
            np = c_loc(mp%H_x)
            dim(1) = mp%g%n1max-mp%g%n1min+1
            dim(2) = mp%nf1max-mp%nf1min+1
#if(_DIM_>=2)        
        case (2)
            np = c_loc(mp%H_y)
            dim(1) = mp%g%n2max-mp%g%n2min+1
            dim(2) = mp%nf2max-mp%nf2min+1
#endif
#if(_DIM_>=3)        
        case (3)
            np = c_loc(mp%H_z)
            dim(1) = mp%g%n3max-mp%g%n3min+1
            dim(2) = mp%nf3max-mp%nf3min+1
#endif
        end select
    end function S(c_get_H_schroedinger)


#endif


#ifndef _REAL_ 
    subroutine S(c_rset_wfs)(psi, f) &
        bind(c, name=SC(rset_wf_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_schroedinger)), pointer :: psip
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
    end subroutine S(c_rset_wfs)
#endif    

    subroutine S(c_set_wfs)(psi, f) &
        bind(c, name=SC(set_wf_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_schroedinger)), pointer :: psip
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
    end subroutine S(c_set_wfs)

#ifndef _REAL_ 
    subroutine S(c_rset_t_wfs)(psi, f, t) &
        bind(c, name=SC(rset_t_wf_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_schroedinger)), pointer :: psip
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
    end subroutine S(c_rset_t_wfs)
#endif    

    subroutine S(c_set_t_wfs)(psi, f, t) &
        bind(c, name=SC(set_t_wf_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_schroedinger)), pointer :: psip
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
    end subroutine S(c_set_t_wfs)



    subroutine S(c_copy_wfs)(psi, source) &
        bind(c, name=SC(copy_wf_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(c_ptr), value :: source 
        type(S(wf_schroedinger)), pointer :: psip, sourcep

        call c_f_pointer(psi, psip)
        call c_f_pointer(source, sourcep)
        call psip%copy(sourcep)
    end subroutine S(c_copy_wfs)

    function S(c_distance_wfs)(psi1, psi2) &
        result(d) bind(c, name=SC(distance_wf_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: psi1 
        type(c_ptr), value :: psi2 
        real(kind = prec) :: d
        type(S(wf_schroedinger)), pointer :: psi1p, psi2p

        call c_f_pointer(psi1, psi1p)
        call c_f_pointer(psi2, psi2p)
        d =  psi1p%distance(psi2p)
    end function S(c_distance_wfs)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! schroedinger specific methods !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine S(c_set_potential_schroedinger)(m, f) &
        bind(c, name=SC(set_potential_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m 
        type(S(schroedinger)), pointer :: mp
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

        call c_f_pointer(m, mp)
        call mp%set_potential(ff)
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
        
    end subroutine S(c_set_potential_schroedinger)

    function S(c_get_potential_schroedinger)(m, dim) &
        result(Vp) bind(c, name=SC(get_potential_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: m 
        type(c_ptr)  :: Vp
        integer, intent(out)  :: dim(_DIM_)
        type(S(schroedinger)), pointer :: mp

        call c_f_pointer(m, mp)
        !Vp = c_loc(mp%V)
#if(_DIM_==1)
        Vp = c_loc(mp%V(mp%g%m1min))
#elif(_DIM_==2)
        Vp = c_loc(mp%V(mp%g%m1min, mp%g%m2min))
#elif(_DIM_==3)
        Vp = c_loc(mp%V(mp%g%m1min, mp%g%m2min, mp%g%m3min))
#endif
        dim(1) = mp%g%m1max-mp%g%m1min+1
#if(_DIM_>=2)        
        dim(2) = mp%g%m2max-mp%g%m2min+1
#endif        
#if(_DIM_>=3)        
        dim(3) = mp%g%m3max-mp%g%m3min+1
#endif        
     end function S(c_get_potential_schroedinger)


    subroutine S(c_load_potential_schroedinger)(m, filename, filename_length) &
        bind(c, name=SC(load_potential_schroedinger)) 
        use iso_c_binding
        type(c_ptr), value :: m 
        integer(c_int), intent(in), value :: filename_length
        character(kind=c_char), dimension(filename_length), intent(in):: filename
        character(len=filename_length, kind=c_char) :: fn 
        type(S(schroedinger)), pointer :: mp
        integer :: j
        call c_f_pointer(m, mp)
        do j=1,filename_length
            fn(j:j) = filename(j)
        end do
        call mp%load_potential(fn)
    end subroutine S(c_load_potential_schroedinger)
 
    subroutine S(c_save_potential_schroedinger)(m, filename, filename_length) &
        bind(c, name=SC(save_potential_schroedinger)) 
        use iso_c_binding
        type(c_ptr), value :: m 
        integer(c_int), intent(in), value :: filename_length
        character(kind=c_char), dimension(filename_length), intent(in):: filename
        character(len=filename_length, kind=c_char) :: fn 
        type(S(schroedinger)), pointer :: mp
        integer :: j
        call c_f_pointer(m, mp)
        do j=1,filename_length
            fn(j:j) = filename(j)
        end do
        call mp%save_potential(fn)
    end subroutine S(c_save_potential_schroedinger)

    subroutine S(c_imaginary_time_propagate_A_wfs)(psi, dt) &
        bind(c, name=SC(imaginary_time_propagate_A_wf_schroedinger)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        type(S(wf_schroedinger)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%imaginary_time_propagate_A(dt)
    end subroutine S(c_imaginary_time_propagate_A_wfs)

    subroutine S(c_imaginary_time_propagate_B_wfs)(psi, dt, method_for_B) &
        bind(c, name=SC(imaginary_time_propagate_B_wf_schroedinger)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        integer, value :: method_for_B
        type(S(wf_schroedinger)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%imaginary_time_propagate_B(dt, method_for_B)
    end subroutine S(c_imaginary_time_propagate_B_wfs)

    subroutine S(c_add_apply_B_wfs)(this, other, coefficient) &
        bind(c, name=SC(add_apply_B_wf_schroedinger)) 
        use iso_c_binding
        type(c_ptr), value :: this
        type(c_ptr), value :: other 
        _COMPLEX_OR_REAL_(kind=prec), value :: coefficient
        type(S(wf_schroedinger)), pointer :: thisp
        type(S(wf_schroedinger)), pointer :: otherp

        call c_f_pointer(this, thisp)
        call c_f_pointer(other, otherp)
        call thisp%add_apply_B(otherp, coefficient)
    end subroutine S(c_add_apply_B_wfs)
   
    function S(c_kinetic_energy_wfs)(psi) &
        result(ans) bind(c, name=SC(kinetic_energy_wf_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_schroedinger)), pointer :: psip
        real(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        ans = psip%kinetic_energy()
    end function S(c_kinetic_energy_wfs)

    function S(c_potential_energy_wfs)(psi) &
        result(ans) bind(c, name=SC(potential_energy_wf_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_schroedinger)), pointer :: psip
        real(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        ans = psip%potential_energy()
    end function S(c_potential_energy_wfs)
    
    function S(c_interaction_energy_wfs)(psi) &
        result(ans) bind(c, name=SC(interaction_energy_wf_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_schroedinger)), pointer :: psip
        real(kind=prec) :: ans

        call c_f_pointer(psi, psip)
        ans = psip%interaction_energy()
    end function S(c_interaction_energy_wfs)

    function S(c_observable_wfs)(psi, f) &
        result(ans) bind(c, name=SC(observable_wf_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: psi 
        type(S(wf_schroedinger)), pointer :: psip
        real(kind=prec) :: ans
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
        ans = psip%observable(ff)
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
    end function S(c_observable_wfs)

    subroutine S(c_get_energy_expectation_deviation_wfs)(psi, ans) &
        bind(c, name=SC(get_energy_expectation_deviation_wf_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: psi 
        real(kind=prec), intent(out) :: ans(2)
        type(S(wf_schroedinger)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%get_energy_expectation_deviation(ans(1), ans(2))
    end subroutine S(c_get_energy_expectation_deviation_wfs)
    
    subroutine S(c_get_realspace_observables_wfs)(psi, ans) &
        bind(c, name=SC(get_realspace_observables_wf_schroedinger))
        use iso_c_binding
        type(c_ptr), value :: psi 
        real(kind=prec), intent(out) :: ans(2+2*_DIM_)
        type(S(wf_schroedinger)), pointer :: psip

        call c_f_pointer(psi, psip)
#if(_DIM_==1)        
        call psip%get_realspace_observables(ans(1), ans(2), ans(3), ans(4))
#elif(_DIM_==2)        
        call psip%get_realspace_observables(ans(1), ans(2), ans(3), ans(4), ans(5), ans(6))
#elif(_DIM_==3)        
        call psip%get_realspace_observables(ans(1), ans(2), ans(3), ans(4), ans(5), ans(6), ans(7), ans(8))
#endif
    end subroutine S(c_get_realspace_observables_wfs)

    subroutine S(c_selfconsistent_nonlinear_step_wfs)(psi, dt, dt1, eps, max_iters) &
        bind(c, name=SC(selfconsistent_nonlinear_step_wf_schroedinger)) 
        use iso_c_binding
        type(c_ptr), value :: psi
        _COMPLEX_OR_REAL_(kind=prec), value :: dt
        _COMPLEX_OR_REAL_(kind=prec), value :: dt1
        real(kind=prec), value :: eps
        integer, value :: max_iters
        type(S(wf_schroedinger)), pointer :: psip

        call c_f_pointer(psi, psip)
        call psip%selfconsistent_nonlinear_step(dt, dt1, eps, max_iters)
    end subroutine S(c_selfconsistent_nonlinear_step_wfs)

#ifdef _QUADPRECISION_
end module S(tssmq_c_schroedinger)
#else
end module S(tssm_c_schroedinger)
#endif
