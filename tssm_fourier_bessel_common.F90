#ifdef _QUADPRECISION_
module tssmq_fourier_bessel_common
use tssmq_base
use tssmq_fourier_common, only: dirichlet, neumann
#else
module tssm_fourier_bessel_common
use tssm_base
use tssm_fourier_common, only: dirichlet, neumann
#endif
implicit none

integer, parameter :: gauss = 1
integer, parameter :: radau = 2
integer, parameter :: lobatto = 3

integer, parameter :: eprec=selected_real_kind(p=30)

contains

function legendre(n, x)
    integer, intent(in) :: n
    real(kind=eprec), intent(in) :: x
    real(kind=eprec) :: legendre
    real(kind=eprec) :: L_k, L_k_minus_1, L_k_plus_1
    integer :: k

    if (n==0) then
        legendre = 1.0_eprec
        return
    elseif (n==1) then
        legendre = x
        return
    endif    
    
    L_k_minus_1 = 1.0_eprec
    L_k = x
    do k = 1, n-1
        L_k_plus_1 = (real(2*k+1, kind=eprec)*x*L_k - real(k, kind=eprec)*L_k_minus_1)/real(k+1, kind=eprec)
        L_k_minus_1 = L_k
        L_k = L_k_plus_1
    end do   
    legendre = L_k
end function legendre



subroutine legendre_deriv(n, x, L_k, Ld_k, forRadau)
    integer, intent(in) :: n
    real(kind=eprec), intent(in) :: x
    real(kind=eprec), intent(out) :: L_k
    real(kind=eprec), intent(out) :: Ld_k
    logical, intent(in), optional :: forRadau

    real(kind=eprec) :: L_k_minus_1, L_k_plus_1
    real(kind=eprec) :: Ld_k_minus_1, Ld_k_plus_1
    integer :: k

    if (n==0) then
        L_k = 1.0_eprec
        Ld_k = 0.0_prec
        return
    elseif (n==1) then
        L_k = x
        Ld_k = 1.0_eprec
        return
    endif    

    L_k_minus_1 = 1.0_eprec
    Ld_k_minus_1 = 0.0_eprec
    L_k = x
    Ld_k = 1.0_eprec 
    do k = 1,n-1
        L_k_plus_1 = (real(2*k+1, kind=eprec)*x*L_k - real(k, kind=eprec)*L_k_minus_1)/real(k+1, kind=eprec)
        Ld_k_plus_1 = (real(2*k+1, kind=eprec)*(x*Ld_k+L_k) - real(k, kind=eprec)*Ld_k_minus_1)/real(k+1, kind=eprec)
        L_k_minus_1 = L_k
        Ld_k_minus_1 = Ld_k
        L_k = L_k_plus_1
        Ld_k = Ld_k_plus_1
    end do   
    if (present(forRadau).and.forRadau) then
       L_k = L_k +  L_k_minus_1 
       Ld_k = Ld_k +  Ld_k_minus_1 
    end if
end subroutine legendre_deriv


subroutine quadrature_coefficients_gauss(x, w, n)
    real(kind=eprec), intent(out) :: x(0:n)
    real(kind=eprec), intent(out) :: w(0:n)
    integer, intent(in) :: n
 
    real(kind=eprec), parameter :: EPS = 10000*epsilon(1.0_eprec)
    real(kind=eprec), parameter :: EPS1 = 10*EPS
    integer, parameter :: MAXIT = 100

    integer :: info
    integer j, k, its
    real(kind=eprec) :: z, L, Lprime

    real(kind=8), allocatable :: d(:), e(:)

    allocate( d(n+1) )
    allocate( e(n) )
    d = 0.0d0
    e = (/ (sqrt(real(j**2, kind=8)/real(4*j**2-1, kind=8)), j = 1,n ) /)

    call DSTERF(n+1, d, e, info)

    if (mod(n,2)==0) then
        x(n/2) = 0.0_eprec
        call legendre_deriv(n+1, 0.0_eprec, L, Lprime)
        w(n/2) = 2.0_eprec/Lprime**2
    end if 

    do k=0, ceiling(0.5*n)-1
        z = real(d(k+1), kind=eprec)
        do its=1, MAXIT
            call legendre_deriv(n+1, z, L, Lprime)
            z = z - L/Lprime
            if (abs(L) <= EPS) then 
                exit
            end if         
        end do
        x(k) = z
        x(n-k) = -z
        call legendre_deriv(n+1, z, L, Lprime)
        w(k) = 2.0_eprec/((1.0_eprec-z**2)*Lprime**2)
        w(n-k) = w(k)
    end do   

    deallocate( d )
    deallocate( e )
    
    ! check that all nodes are distinct and increasing
    do k=1,n
        if(.not.((x(k)-x(k-1))>EPS1)) then
            stop 'nodes not distinct in gauss_lobatto'
        end if
    end do    
end subroutine quadrature_coefficients_gauss


subroutine quadrature_coefficients_gauss_radau(x, w, n)
    real(kind=eprec), intent(out) :: x(0:n)
    real(kind=eprec), intent(out) :: w(0:n)
    integer, intent(in) :: n
 
    real(kind=eprec), parameter :: EPS = 10000*epsilon(1.0_eprec)
    real(kind=eprec), parameter :: EPS1 = 10*EPS
    integer, parameter :: MAXIT = 100

    integer :: info
    integer j, k, its
    real(kind=eprec) :: z, L, Lprime

    real(kind=8), allocatable :: d(:), e(:)

    allocate( d(n) )
    allocate( e(n-1) )
    d = (/ (1.0_prec/real((2*j+1)*(2*j+3), kind=8), j=0,n-1 ) /)
    e = (/ (sqrt(real(j*(j+1), kind=8)/real((2*j+1)**2, kind=8)), j = 1,n-1 ) /)

    call DSTERF(n, d, e, info)

    x(0) = -1.0_eprec
    w(0) = 2.0_eprec/(real((n+1)**2, kind=eprec)*legendre(n+1, -1.0_eprec)**2)
    do k=1,n 
        z = real(d(k), kind=eprec)
        do its=1, MAXIT
            call legendre_deriv(n+1, z, L, Lprime, forRadau=.true.)
            z = z - L/Lprime
            if (abs(L) <= EPS) then 
                exit
            end if         
        end do
        x(k) = z
        w(k) = (1.0_eprec-z)/(real((n+1)**2, kind=eprec)*legendre(n+1, z)**2)
    end do   

    deallocate( d )
    deallocate( e )
    
    ! check that all nodes are distinct and increasing
    do k=1,n
        if(.not.((x(k)-x(k-1))>EPS1)) then
            stop 'nodes not distinct in gauss_radau'
        end if
    end do    
end subroutine quadrature_coefficients_gauss_radau



subroutine quadrature_coefficients_gauss_lobatto(x, w, n)
    real(kind=eprec), intent(out) :: x(0:n)
    real(kind=eprec), intent(out) :: w(0:n)
    integer, intent(in) :: n
 
    real(kind=eprec), parameter :: EPS = 10000*epsilon(1.0_eprec)
    real(kind=eprec), parameter :: EPS1 = 10*EPS
    integer, parameter :: MAXIT = 100

    integer :: info
    integer j, k, its
    real(kind=eprec) :: z, L, Lprime

    real(kind=8), allocatable :: d(:), e(:)

    allocate( d(n-1) )
    allocate( e(n-2) )
    d = 0.0d0
    e = (/ (sqrt(real(j*(j+2), kind=8)/real((2*j+1)*(2*j+3), kind=8)), j = 1,n-2 ) /)

    call DSTERF(n-1, d, e, info)

    x(0) = -1.0_eprec
    x(n) = 1.0_eprec
    w(n) = 2.0_eprec/(real(n*(n+1), kind=eprec)*legendre(n, 1.0_eprec)**2)
    w(0) = w(n)
    if (mod(n,2)==0) then
        x(n/2) = 0.0_eprec
        w(n/2) = 2.0_eprec/(real(n*(n+1), kind=eprec)*legendre(n, 0.0_eprec)**2)
    end if 

    do k=1, ceiling(0.5*n)-1
        z = real(d(k), kind=eprec)
        do its=1, MAXIT
            call legendre_deriv(n, z, L, Lprime)
            z = z - (1.0_eprec - z**2)*Lprime/(2.0_eprec*z*Lprime - real(n*(n+1),kind=eprec)*L)
            if (abs(Lprime) <= EPS) then 
                exit
            end if         
        end do
        x(k) = z
        x(n-k) = -z
        w(k) = 2.0_eprec/(real(n*(n+1), kind=eprec)*legendre(n, z)**2)
        w(n-k) = w(k)
    end do   

    deallocate( d )
    deallocate( e )

    ! check that all nodes are distinct and increasing
    do k=1,n
        if(.not.((x(k)-x(k-1))>EPS1)) then
            stop 'nodes not distinct in gauss_lobatto'
        end if
    end do    
end subroutine quadrature_coefficients_gauss_lobatto


function besselj_zero_approx1(nu, m)
    integer, intent(in) :: nu
    integer, intent(in) :: m
    real(kind=prec) :: besselj_zero_approx1

    real(kind=prec), parameter :: pi = 3.1415926535897932384626433832795028841971693993751_prec
    real(kind=prec) :: mu, a, a8, a82
    
    mu = 4.0_prec*real(nu, kind=prec)**2
    a = (real(m, kind=prec)+0.5_prec*real(nu, kind=prec)-0.25_prec)*pi
    a8 = 1.0_prec/(8.0_prec*a)
    a82 = a8**2
    besselj_zero_approx1  = a - (mu-1.0_prec)*a8*(1.0_prec + a82*(4.0_prec*(7.0_prec*mu-31.0_prec)/3.0_prec &
         +a82*( 32.0_prec*(83.0_prec*mu**2-982.0_prec*mu+3779.0_prec)/15.0_prec &
         +a82*(64.0_prec*(6949.0_prec*mu**3-153855.0_prec*mu**2+1585743.0_prec*mu-6277237.0_prec)/105.0_prec))))
end function besselj_zero_approx1


function besseljprime_zero_approx1(nu, m)
    integer, intent(in) :: nu
    integer, intent(in) :: m
    real(kind=prec) :: besseljprime_zero_approx1

    real(kind=prec), parameter :: pi = 3.1415926535897932384626433832795028841971693993751_prec
    real(kind=prec) :: m1, mu, b, b8, b82
    if (nu==0) then
        m1 = m +1
    else
        m1 = m
    endif        
    mu = 4.0_prec*real(nu, kind=prec)**2
    b = (real(m1, kind=prec)+0.5_prec*real(nu, kind=prec)-0.75_prec)*pi
    b8 = 1.0_prec/(8.0_prec*b)
    b82 = b8**2
    besseljprime_zero_approx1 = & 
           b - b8*(mu+3 + b82*(4.0_prec*(7.0_prec*mu**2+82.0_prec*mu-9)/3.0_prec &
           +b82*(32.0_prec*(83.0_prec*mu**3+2075.0_prec*mu**2-3039.0_prec*mu+3537.0_prec)/15.0_prec &
           +b82*(64.0_prec*(6949.0_prec*mu**4+296492.0_prec*mu**3-1248002.0_prec*mu**2+7414380.0_prec*mu-585362)/105.0_prec))))
end function besseljprime_zero_approx1


function airyAi_zero_approx(m)
    integer, intent(in) :: m
    real(kind=prec) :: airyAi_zero_approx
    real(kind=prec) :: t, t2
    real(kind=prec), parameter :: pi = 3.1415926535897932384626433832795028841971693993751_prec
    real(kind=prec), dimension(10), parameter :: small_airy_zeros = &
       (/ -2.338107410459767_prec, -4.087949444130970_prec, -5.520559828095552_prec, &
          -6.786708090071759_prec, -7.944133587120853_prec, -9.022650853340980_prec, &
          -10.04017434155809_prec, -11.00852430373326_prec, -11.93601556323626_prec, &
          -12.82877675286576_prec /)

    if (m<=10) then
        airyAi_zero_approx = small_airy_zeros(m)
    else
        t = 0.375_prec*pi*(4.0_prec*m-1.0_prec)
        t2 = t**(-2)
        airyAi_zero_approx = -t**(2.0_prec/3.0_prec)*(1.0_prec+t2*(5.0_prec/48.0_prec &
            +t2*(-5.0_prec/36.0_prec+t2*(77125.0_prec/82944.0_prec &
            +t2*(-108056875.0_prec/6967296.0_prec+t2*162375596875.0_prec/334430208.0_prec))))) 
    endif
end function airyAi_zero_approx


function airyAiprime_zero_approx(m)
    integer, intent(in) :: m
    real(kind=prec) :: airyAiprime_zero_approx
    real(kind=prec) :: t, t2
    real(kind=prec), parameter :: pi = 3.1415926535897932384626433832795028841971693993751_prec
    real(kind=prec), dimension(10), parameter :: small_airyprime_zeros = &
      (/ -1.0187929716474710890, -3.2481975821798365379, -4.8200992111787356394, &
         -6.1633073556394865476, -7.3721772550477701771, -8.4884867340197221329, &
         -9.5354490524335474707, -10.527660396957407282, -11.475056633480245295, &
         -12.384788371845747325 /)
    if (m<=10) then
        airyAiprime_zero_approx = small_airyprime_zeros(m)
    else
        t = 0.375_prec*pi*(4.0_prec*m-3.0_prec)
        t2 = t**(-2)
        airyAiprime_zero_approx = -t**(2.0_prec/3.0_prec)*(1.0_prec+t2*(-7.0_prec/48.0_prec &
            +t2*(35.0_prec/288.0_prec+t2*(-181223.0_prec/207360.0_prec &
            +t2*(18683371.0_prec/1244160.0_prec-t2*9114588436.0_prec/191102976.0_prec))))) 
    endif
end function airyAiprime_zero_approx


function besselj_zero_approx2(nu, m, prime)
    integer, intent(in) :: nu
    integer, intent(in) :: m
    logical, intent(in), optional :: prime
    real(kind=prec) :: besselj_zero_approx2
    real(kind=prec), parameter :: pi = 3.1415926535897932384626433832795028841971693993751_prec
    real(kind=prec) :: z_airy, zeta, y, x, p, p2, x2, r, z, h, f1, g1

    if (present(prime).and.prime) then
        z_airy = airyAiprime_zero_approx(m)
    else
        z_airy = airyAi_zero_approx(m)
    end if    
    zeta = nu**(-2.0_prec/3.0_prec)*z_airy
    y = 2.0_prec/3.0_prec*(-zeta)**(3.0_prec/2.0_prec)
    if (y>100000.0_prec) then
        x = 0.5_prec*pi
    elseif (y<1.0_prec)then
        p = (3.0_prec*y)**(1.0_prec/3.0_prec)
        p2 = p**2
        x = p*(1.0_prec+p2*(-2.0_prec/15.0_prec+p2*(3.0_prec/175.0_prec+p2*(-2.0_prec/1575.0_prec))))
    else
        p = 1.0_prec/(y+0.5_prec*pi)
        p2 = p**2
        x = 0.5*pi - p*(1.0_prec+p2*(2.0_prec/3.0_prec+p2*(13.0_prec/15.0_prec &
             +p2*(146.0_prec/105.0_prec+p2*(781.0_prec/315.0_prec+p2*16328.0_prec/3465.0_prec)))))
    endif
    x2 = (y+x)**2
    r = (x-atan(x+y))/x2
    x = x - (1.0_prec+x2)*r*(1.0_prec+r/(x+y))

    z = 1.0_prec/cos(x)
    h = sqrt(zeta*(1.0_prec-z**2))
    if (present(prime).and.prime) then
        g1 = -z/zeta*h*(7.0_prec/(48.0_prec*zeta) + h*(7.0_prec/(z**2-1.0_prec)+9.0_prec)/24.0_prec)
        besselj_zero_approx2 = nu*z + g1/nu
    else
        besselj_zero_approx2 = nu*z + g1/nu
        f1 = -z/zeta*h*(5.0_prec/(48.0_prec*zeta) + h*(5.0_prec/(z**2-1.0_prec)+3.0_prec)/24.0_prec)
        besselj_zero_approx2 = nu*z + f1/nu
    end if
end function besselj_zero_approx2


function besselj_zero_approx(nu, m)
    integer, intent(in) :: nu
    integer, intent(in) :: m
    real(kind=prec) :: besselj_zero_approx

    if (m>=nu) then
        besselj_zero_approx = besselj_zero_approx1(nu, m)
    else
        besselj_zero_approx = besselj_zero_approx2(nu, m)
    endif
end function besselj_zero_approx


function besseljprime_zero_approx(nu, m)
    integer, intent(in) :: nu
    integer, intent(in) :: m
    real(kind=prec) :: besseljprime_zero_approx

    if (m>=nu) then
        besseljprime_zero_approx = besseljprime_zero_approx1(nu, m)
    else
        besseljprime_zero_approx = besselj_zero_approx2(nu, m, prime=.true.)
    endif
end function besseljprime_zero_approx


function besselj_zero_iter(nu, z0) result(z)
    integer, intent(in) :: nu
    real(kind=eprec), intent(in) :: z0
    real(kind=eprec) :: z

    real(kind=eprec), parameter :: EPS = 10000*epsilon(1.0_eprec)
    integer, parameter :: MAXITER = 100

    integer :: its
    real(kind=eprec) :: J, Jprime

    z = z0

    do its = 1, MAXITER
        J = bessel_jn(nu, z)
        Jprime = bessel_jn(nu-1, z) - real(nu, kind=eprec)*J/z
        z = z - J/Jprime   
        if (abs(J) <= EPS) then
            exit
        end if       
    end do
end function besselj_zero_iter


function besselj_zero(nu, m)
    integer, intent(in) :: nu
    integer, intent(in) :: m
    real(kind=eprec) :: besselj_zero

    besselj_zero = besselj_zero_iter(nu, real(besselj_zero_approx(nu, m), kind=eprec))
end function besselj_zero


function besseljprime_zero_iter(nu, z0) result(z)
    integer, intent(in) :: nu
    real(kind=eprec), intent(in) :: z0
    real(kind=eprec) :: z

    real(kind=eprec), parameter :: EPS = 10000*epsilon(1.0_eprec)
    integer, parameter :: MAXITER = 100

    integer :: its
    real(kind=eprec) :: J, J1, J2

    z = z0

    do its = 1, MAXITER
        J = bessel_jn(nu, z) 
        J1 = bessel_jn(nu-1, z) - real(nu, kind=eprec)*J/z
        J2 = -J1/z + ((nu/z)**2 - 1.0_eprec)*J
        z = z - J1/J2
        if (abs(J1) <= EPS) then
            exit
        end if       
    end do
end function besseljprime_zero_iter


function besseljprime_zero(nu, m)
    integer, intent(in) :: nu
    integer, intent(in) :: m
    real(kind=eprec) :: besseljprime_zero

    besseljprime_zero = besseljprime_zero_iter(nu, real(besseljprime_zero_approx(nu, m), kind=eprec))
end function besseljprime_zero



#ifdef _QUADPRECISION_
subroutine fourier_bessel_coeffs(nr, kk, mm, x, w, L, eigenvalues, normalization_factors, nfthetamin, nfthetamax, &
                                 boundary_conditions, quadrature_formula)
#else
subroutine fourier_bessel_coeffs_eprec(nr, kk, mm, x, w, L, eigenvalues, normalization_factors, nfthetamin, nfthetamax, &
                                       boundary_conditions, quadrature_formula)
#endif
    integer, intent(in) :: nr, kk, mm
    real(kind=eprec), intent(out) :: x(1:nr)
    real(kind=eprec), intent(out) :: w(1:nr)
    real(kind=eprec), intent(out) :: L(1:nr, 1:kk, 0:mm/2) 
    real(kind=eprec), intent(out) :: eigenvalues(1:kk, nfthetamin:nfthetamax)
    real(kind=eprec), intent(out) :: normalization_factors(1:kk, nfthetamin:nfthetamax)
    integer, intent(in) :: nfthetamin, nfthetamax 
    integer, intent(in) :: boundary_conditions
    integer, intent(in) :: quadrature_formula

    integer :: m, k, n
    real(kind=eprec) :: z
    real(kind=eprec) :: z2
    real(kind=eprec) :: f 
    real(kind=eprec), parameter :: pi = 3.1415926535897932384626433832795028841971693993751_eprec
    
    real(kind=eprec), allocatable :: xx(:), ww(:)

    allocate( xx(0:nr+1) )
    allocate( ww(0:nr+1) )

    select case(quadrature_formula)
    case (gauss) 
        call quadrature_coefficients_gauss(xx, ww, nr-1)
        x = 0.5_eprec*xx(0:nr-1) + 0.5_eprec ! transform [-1,1] -> [0,1]
        w = 0.5_eprec*ww(0:nr-1)*x
    case (radau) 
        call quadrature_coefficients_gauss_radau(xx, ww, nr)
        x = 0.5_eprec*xx(1:nr) + 0.5_eprec ! transform [-1,1] -> [0,1]
        w = 0.5_eprec*ww(1:nr)*x
    case (lobatto) 
        select case(boundary_conditions)
        case (dirichlet)
            call quadrature_coefficients_gauss_lobatto(xx, ww, nr+1)
        case (neumann)
            call quadrature_coefficients_gauss_lobatto(xx, ww, nr)
        end select    
        x = 0.5_eprec*xx(1:nr) + 0.5_eprec ! transform [-1,1] -> [0,1]
        w = 0.5_eprec*ww(1:nr)*x
    end select

    deallocate( xx, ww )
!$OMP PARALLEL DO PRIVATE(m, k, z, z2, f)
    do m=0,mm/2
        do k=1,kk
           select case(boundary_conditions)
           case (dirichlet)
               z = besselj_zero(m, k)
               z2 = z**2
               f = sqrt(1.0_eprec/pi)/abs(bessel_jn(m+1, z)) 
           case (neumann)
               z = besseljprime_zero(m, k)
               z2 = z**2
               f = sqrt(1.0_eprec/(pi*(1.0_eprec - real(m**2,kind=eprec)/z2)))/abs(bessel_jn(m, z)) 
           end select   
           if ((m>=nfthetamin).and.(m<=nfthetamax)) then
               eigenvalues(k, m) = -z2
               normalization_factors(k, m) = f
           end if
           if ((m>=1).and.(m<mm/2).and.(mm-m>=nfthetamin).and.(mm-m<=nfthetamax)) then
               eigenvalues(k, mm-m) = -z2
               normalization_factors(k, mm-m) = f
           end if           
           do n = 1,nr
               L(n,k,m) = f*bessel_jn(m, z*x(n))
           end do
        end do
    end do    
!$OMP END PARALLEL DO

#ifdef _QUADPRECISION_
end subroutine fourier_bessel_coeffs
#else
end subroutine fourier_bessel_coeffs_eprec
#endif


#ifndef _QUADPRECISION_
subroutine fourier_bessel_coeffs(nr, kk, mm, x, w, L, eigenvalues, normalization_factors, nfthetamin, nfthetamax, &
                                 boundary_conditions, quadrature_formula)
    integer, intent(in) :: nr, kk, mm
    real(kind=prec), intent(out) :: x(1:nr)
    real(kind=prec), intent(out) :: w(1:nr)
    real(kind=prec), intent(out) :: L(1:nr, 1:kk, 0:mm/2) 
    real(kind=prec), intent(out) :: eigenvalues(1:kk, nfthetamin:nfthetamax)
    real(kind=prec), intent(out) :: normalization_factors(1:kk, nfthetamin:nfthetamax)
    integer, intent(in) :: nfthetamin, nfthetamax 
    integer, intent(in) :: boundary_conditions
    integer, intent(in) :: quadrature_formula

    real(kind=eprec) :: x_eprec(1:nr)
    real(kind=eprec) :: w_eprec(1:nr)
    real(eprec), allocatable :: L_eprec(:,:,:)
    real(eprec), allocatable :: ev_eprec(:,:)
    real(eprec), allocatable :: nf_eprec(:,:)

    allocate ( L_eprec(1:nr, 1:kk, 0:mm/2) ) 
    allocate ( ev_eprec(1:kk, nfthetamin:nfthetamax) ) 
    allocate ( nf_eprec(1:kk, nfthetamin:nfthetamax) ) 

    call fourier_bessel_coeffs_eprec(nr, kk, mm, x_eprec, w_eprec, L_eprec, &
                                     ev_eprec, nf_eprec, nfthetamin, nfthetamax, &
                                     boundary_conditions, quadrature_formula)

    x = real(x_eprec, kind=prec) ! conversion eprec->double should be done implicitely!!!
    w = real(w_eprec, kind=prec)
    L = real(L_eprec, kind=prec)
    eigenvalues = real(ev_eprec, kind=prec)
    normalization_factors = real(nf_eprec, kind=prec)

    deallocate (L_eprec)
    deallocate (ev_eprec)
    deallocate (nf_eprec)
end subroutine fourier_bessel_coeffs
#endif


#ifdef _QUADPRECISION_
subroutine bessel_rotsym_coeffs(nr, x, w, L, eigenvalues, normalization_factors, boundary_conditions)
#else
subroutine bessel_rotsym_coeffs_eprec(nr, x, w, L, eigenvalues, normalization_factors, boundary_conditions)
#endif
    integer, intent(in) :: nr
    real(kind=eprec), intent(out) :: x(1:nr)
    real(kind=eprec), intent(out) :: w(1:nr)
    real(kind=eprec), intent(out) :: L(1:nr, 1:nr) 
    real(kind=eprec), intent(out) :: eigenvalues(1:nr)
    real(kind=eprec), intent(out) :: normalization_factors(1:nr)
    integer, intent(in) :: boundary_conditions

    real(kind=eprec), parameter :: pi = 3.1415926535897932384626433832795028841971693993751_eprec
    integer :: k, m
    real(kind=eprec) :: z_max, f
    real(kind=eprec), allocatable :: z(:)

    allocate( z(1:nr) )

!TODO: Neumann boundary conditions!

    z_max = besselj_zero(0, nr+1) 
    do k=1,nr
        z(k) = besselj_zero(0, k)
        x(k) = z(k)/z_max
        w(k) = 4.0_eprec*pi/(z_max*bessel_j1(z(k)))**2
        eigenvalues(k) = -z(k)**2
    end do
    do k=1,nr
        f = sqrt(1.0_eprec/pi)/abs(bessel_j1(z(k))) 
        normalization_factors(k) = f
        do m = 1,nr
            L(m,k) = f*bessel_j0(z(k)*x(m))
        end do
    end do

    deallocate( z )

#ifdef _QUADPRECISION_
end subroutine bessel_rotsym_coeffs
#else
end subroutine bessel_rotsym_coeffs_eprec
#endif


#ifndef _QUADPRECISION_
subroutine bessel_rotsym_coeffs(nr, x, w, L, eigenvalues, normalization_factors, boundary_conditions)
    integer, intent(in) :: nr
    real(kind=prec), intent(out) :: x(1:nr)
    real(kind=prec), intent(out) :: w(1:nr)
    real(kind=prec), intent(out) :: L(1:nr, 1:nr) 
    real(kind=prec), intent(out) :: eigenvalues(1:nr)
    real(kind=prec), intent(out) :: normalization_factors(1:nr)
    integer, intent(in) :: boundary_conditions

    real(eprec), allocatable :: x_eprec(:)
    real(eprec), allocatable :: w_eprec(:)
    real(eprec), allocatable :: L_eprec(:,:)
    real(eprec), allocatable :: ev_eprec(:)
    real(eprec), allocatable :: nf_eprec(:)

    allocate ( x_eprec(1:nr) ) 
    allocate ( w_eprec(1:nr) ) 
    allocate ( L_eprec(1:nr, 1:nr) ) 
    allocate ( ev_eprec(1:nr) ) 
    allocate ( nf_eprec(1:nr) ) 

    call bessel_rotsym_coeffs_eprec(nr, x_eprec, w_eprec, L_eprec, ev_eprec, &
                                    nf_eprec, boundary_conditions)

    x = real(x_eprec, kind=prec) 
    w = real(w_eprec, kind=prec)
    L = real(L_eprec, kind=prec)
    eigenvalues = real(ev_eprec, kind=prec)
    normalization_factors = real(nf_eprec, kind=prec)

    deallocate (x_eprec)
    deallocate (w_eprec)
    deallocate (L_eprec)
    deallocate (ev_eprec)
end subroutine bessel_rotsym_coeffs
#endif


#ifdef _QUADPRECISION_
end module tssmq_fourier_bessel_common
#else
end module tssm_fourier_bessel_common
#endif
