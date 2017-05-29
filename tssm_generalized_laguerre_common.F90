#ifdef _QUADPRECISION_
module tssmq_generalized_laguerre_common
use tssmq_base
#else
module tssm_generalized_laguerre_common
use tssm_base
#endif
implicit none

integer, parameter :: eprec=selected_real_kind(p=30)

contains

subroutine  gauss_laguerre(x, w, n, alpha)
    real(kind=eprec), intent(out) :: x(n)
    real(kind=eprec), intent(out) :: w(n)
    integer, intent(in) :: n
    real(kind=eprec), intent(in), optional :: alpha

    real(kind=eprec), parameter :: EPS = 10000*epsilon(1.0_eprec)
    real(kind=eprec), parameter :: EPS1 = 10*EPS
    integer, parameter :: MAXIT = 30

    real(kind=eprec) :: z, z1, p1, p2, p3, pp, alf
    integer :: i, j, its, ai

    if (present(alpha)) then
       alf = alpha
    else
       alf = 0.0_eprec
    end if   

    do i=1,n
    if (i==1) then
        z = (1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf)
    elseif (i==2) then
        z = z + (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n)
    else 
        ai = i-2;
        z = z + ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf / &
                (1.0+3.5*ai))*(z-x(i-2))/(1.0+0.3*alf)
    end if
    do its=1,MAXIT
        p1 = 1.0
        p2 = 0.0
        do j=1,n
            p3 = p2
            p2 = p1
            p1 = ((real(2*j-1, kind=eprec)+alf-z)*p2-(real(j-1, kind=eprec)+alf)*p3)/real(j, kind=eprec)
        end do
        pp = (real(n, kind=eprec)*p1-(real(n, kind=eprec)+alf)*p2)/z
        z1 = z
        z = z1-p1/pp;
        if (abs(z-z1) <= EPS) then
            exit
        end if
    end do
    if (its > MAXIT) then
        print *, 'delta =', abs(z-z1)
        print *, 'eps =', EPS
        print *, 'maxit =', MAXIT 
        stop 'too many iterations in gauss_laguerre'
    end if
    x(i)=z
    w(i) = -exp(log_gamma(alf+n)-log_gamma(real(n, kind=eprec)))/(pp*n*p2)
    end do

    ! check that all nodes are distinct and increasing
    do i=2,n
        if(.not.((x(i)-x(i-1))>EPS1)) then
            stop 'nodes not distinct in gauss_laguerre'
        end if
    end do    

end subroutine  gauss_laguerre

subroutine gauss_laguerre_scaled(x, w, n, gg, m)
    real(kind=eprec), intent(out) :: x(n)
    real(kind=eprec), intent(out) :: w(n)
    integer, intent(in) :: n
    real(kind=eprec), intent(in) :: gg
    integer, intent(in), optional :: m

    real(kind=eprec) :: one_over_gg
    integer :: i

    one_over_gg = 1.0_eprec/gg

    if (present(m).and.m/=0) then
        call gauss_laguerre(x, w, n, real(m, kind=eprec))
        do i=1,n
            w(i) = w(i)*exp(x(i))*0.5_eprec*one_over_gg*x(i)**(-m)
        end do
    else
        call gauss_laguerre(x, w, n)
        do i=1,n
            w(i) = w(i)*exp(x(i))*0.5_eprec*one_over_gg
        end do
    end if

    do i=1,n
        x(i) = sqrt(x(i)*one_over_gg)
    end do

end subroutine gauss_laguerre_scaled

subroutine laguerre_polynomial(x, n, L, alpha)
    real(kind=eprec), intent(in) :: x 
    integer, intent(in) :: n
    real(kind=eprec), intent(out) :: L(0:n)
    real(kind=eprec), intent(in), optional :: alpha

    real(kind=eprec) :: alpha_minus_x 
    integer :: i

    if (n>=0) then
        L(0) = 1.0_eprec
    end if
    if (n>=1) then
       if (present(alpha)) then
           alpha_minus_x = alpha - x
           L(1) = alpha - x + 1.0_eprec
           do i=2,n
               L(i) = ((alpha_minus_x + real(2*i-1, kind=eprec))*L(i-1) - (alpha+real(i-1, kind=eprec))*L(i-2)) &
            &         /real(i, kind=eprec)
           end do
       else
           L(1) = 1.0_eprec - x
           do i=2,n
               L(i) = ((real(2*i-1, kind=eprec)-x)*L(i-1) - real(i-1, kind=eprec)*L(i-2)) &
            &         /real(i, kind=eprec)
           end do
       end if
    end if
end subroutine laguerre_polynomial



#ifdef _QUADPRECISION_
subroutine laguerre_scaled_coeffs(n, x, w, L, gg, m)
#else
subroutine laguerre_scaled_coeffs_eprec(n, x, w, L, gg, m)
#endif
    integer, intent(in) :: n
    real(kind=eprec), intent(out) :: x(0:n)
    real(kind=eprec), intent(out) :: w(0:n)
    real(kind=eprec), intent(out) :: L(0:n, 0:n)
    real(kind=eprec), intent(in) :: gg 
    integer, intent(in), optional :: m
    
    real(kind=eprec), parameter :: pi = 3.1415926535897932384626433832795028841971693993751_eprec
    real(kind=eprec) :: sqrt_gg_over_pi, sqrt_gg
    integer :: k, j 
    real(kind=eprec) :: ll(0:n)
    real(kind=eprec) :: f_j 

    real(kind=eprec) :: cc(0:n)

    sqrt_gg_over_pi = sqrt(gg/pi)

    if (present(m).and.m/=0) then
        call gauss_laguerre_scaled(x, w, n+1, gg, m)

        cc(0) = gamma(real(m+1, kind=eprec)) 
        do k=1,n
            cc(k) = cc(k-1)*real(m+k, kind=eprec)/real(k, kind=eprec)
        end do  
        do k=0,n
            cc(k) = 1.0_eprec/sqrt(cc(k))
        end do  
        sqrt_gg = sqrt(gg)

        do j=0,n
            f_j = exp(-0.5_eprec*gg*x(j)**2)*sqrt_gg_over_pi*(sqrt_gg*x(j))**m
            call laguerre_polynomial(gg*x(j)**2, n, ll, real(m, kind=eprec))
            do k=0,n
                L(j,k) = f_j*cc(k)*ll(k)
            end do
        end do
    else
        call gauss_laguerre_scaled(x, w, n+1, gg)

        do j=0,n
            f_j = exp(-0.5_eprec*gg*x(j)**2)*sqrt_gg_over_pi
            call laguerre_polynomial(gg*x(j)**2, n, ll)
            do k=0,n
                L(j,k) = f_j*ll(k)
            end do
        end do    
    end if
#ifdef _QUADPRECISION_
end subroutine laguerre_scaled_coeffs
#else
end subroutine laguerre_scaled_coeffs_eprec
#endif

#ifndef _QUADPRECISION_
subroutine laguerre_scaled_coeffs_prec(n, x, w, L, gg, m)
    integer, intent(in) :: n
    real(kind=prec), intent(out) :: x(0:n)
    real(kind=prec), intent(out) :: w(0:n)
    real(kind=prec), intent(out) :: L(0:n, 0:n)
    real(kind=prec), intent(in) :: gg 
    integer, intent(in), optional :: m

    real(eprec) :: x_eprec(0:n)
    real(eprec) :: w_eprec(0:n)
    real(eprec), allocatable :: L_eprec(:,:)

    allocate (L_eprec(0:n, 0:n))

    if (present(m).and.m/=0) then
        call laguerre_scaled_coeffs_eprec(n, x_eprec, w_eprec, L_eprec, real(gg, kind=eprec), m)
    else    
        call laguerre_scaled_coeffs_eprec(n, x_eprec, w_eprec, L_eprec, real(gg, kind=eprec))
    end if    

    x = real(x_eprec, kind=prec) ! conversion eprec->c_double should be done implicitely!!!
    w = real(w_eprec, kind=prec)
    L = real(L_eprec, kind=prec)

    deallocate (L_eprec)
end subroutine laguerre_scaled_coeffs_prec
#endif



#ifdef _QUADPRECISION_
subroutine generalized_laguerre_scaled_coeffs(kk, mm, x, w, L, gg)
#else
subroutine generalized_laguerre_scaled_coeffs_eprec(kk, mm, x, w, L, gg)
#endif
    integer, intent(in) :: kk, mm
    real(kind=eprec), intent(out) :: x(0:kk+mm/2)
    real(kind=eprec), intent(out) :: w(0:kk+mm/2)
    real(kind=eprec), intent(out) :: L(0:kk+mm/2, 0:kk, 0:mm/2) 
    real(kind=eprec), intent(in) :: gg 
    
    real(kind=eprec), parameter :: pi = 3.1415926535897932384626433832795028841971693993751_eprec
    real(kind=eprec) :: sqrt_gg, sqrt_gg_over_pi 
    integer :: k, j, m 
    real(kind=eprec) :: ll(0:kk) 
    real(kind=eprec) :: cc(0:kk, 0:mm/2) 
    real(kind=eprec) :: f_j, f_m

    do k=0,kk 
        cc(k,0) = 1.0_eprec 
        do m=1,mm/2
            cc(k,m) = cc(k,m-1) * real(k+m, kind=eprec)
        end do  
    end do    
    do k=0,kk 
        do m=0,mm/2
            cc(k,m) = 1.0_eprec/sqrt(cc(k,m))
        end do
    end do    

    sqrt_gg = sqrt(gg)
    sqrt_gg_over_pi = sqrt(gg/pi)

    call gauss_laguerre_scaled(x, w, kk+mm/2+1, gg)

    do j=0,kk+mm/2
        f_j = exp(-0.5_eprec*gg*x(j)**2)*sqrt_gg_over_pi
        do m=0,mm/2
            f_m = f_j*(sqrt_gg*x(j))**m
            call laguerre_polynomial(gg*x(j)**2, kk, ll, real(m, kind=eprec))  
            do k=0,kk 
                L(j,k,m) = f_m*cc(k,m)*ll(k)
            end do
        end do    
    end do    
#ifdef _QUADPRECISION_
end subroutine generalized_laguerre_scaled_coeffs
#else
end subroutine generalized_laguerre_scaled_coeffs_eprec
#endif


#ifndef _QUADPRECISION_
subroutine generalized_laguerre_scaled_coeffs(kk, mm, x, w, L, gg)
    integer, intent(in) :: kk, mm
    real(kind=prec), intent(out) :: x(0:kk+mm/2)
    real(kind=prec), intent(out) :: w(0:kk+mm/2)
    real(kind=prec), intent(out) :: L(0:kk+mm/2, 0:kk, 0:mm/2) 
    real(kind=prec), intent(in) :: gg 

    real(kind=eprec) :: x_eprec(0:kk+mm/2)
    real(kind=eprec) :: w_eprec(0:kk+mm/2)
    real(eprec), allocatable :: L_eprec(:,:,:)

    allocate (L_eprec(0:kk+mm/2, 0:kk, 0:mm/2)) 

    call generalized_laguerre_scaled_coeffs_eprec(kk, mm, x_eprec, w_eprec, L_eprec, real(gg, kind=eprec))

    x = real(x_eprec, kind=prec) ! conversion eprec->double should be done implicitely!!!
    w = real(w_eprec, kind=prec)
    L = real(L_eprec, kind=prec)

    deallocate (L_eprec)
end subroutine generalized_laguerre_scaled_coeffs
#endif


#ifdef _QUADPRECISION_
end module tssmq_generalized_laguerre_common
#else
end module tssm_generalized_laguerre_common
#endif
