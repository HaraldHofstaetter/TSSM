module tssm_hermite_common
use tssm
implicit none

integer, parameter :: eprec=selected_real_kind(p=30)

contains

subroutine gauss_hermite(x, w, n)
    real(kind=eprec), intent(out) :: x(n)
    real(kind=eprec), intent(out) :: w(n)
    integer, intent(in) :: n
 
    real(kind=eprec), parameter :: EPS = 10000*epsilon(1.0_eprec)
    real(kind=eprec), parameter :: EPS1 = 10*EPS
    integer, parameter :: MAXIT = 100
    real(kind=eprec), parameter :: PIM4 = 0.7511255444649424828587030047762276930524_eprec ! 1/pi**(1/4)

    real(kind=eprec) :: z, z1, p1, p2, p3, pp
    integer :: m, i, j, its

    real(kind=8), allocatable :: d(:), e(:)
    integer :: info

    allocate( d(n) )
    allocate( e(n-1) )
    d = 0.0d0
    e = (/ (sqrt(0.5d0*j), j=1,n-1) /)

    call DSTERF(n, d, e, info)

    m = (n+1)/2

    do i=1,m
    z = real(d(i), kind=eprec)
!    if (i==1) then
!        z = sqrt(real(2*n+1))-1.85575*(2*n+1)**(-0.16667)
!    elseif (i==2) then
!        z = -x(1) + 1.14*(n**0.426)/x(1)
!    elseif (i==3) then
!        z = -1.86*x(2)+0.86*x(1)
!    elseif (i==4) then
!        z = -1.91*x(3)+0.91*x(2)
!    else
!        z = -2.0*x(i-1)+x(i-2)
!    end if
    do its=1,MAXIT
        p1=PIM4
        p2=0.0
        do j=1,n
            p3=p2
            p2=p1
            p1=z*sqrt(2.0/real(j, kind=eprec))*p2-sqrt(((j-1))/real(j, kind=eprec))*p3
        end do
        pp = sqrt(real(2*n, kind=eprec))*p2
        z1 = z
        z = z1-p1/pp
        if (abs(z-z1) <= EPS) then
            exit 
        end if
    end do
    if (its > MAXIT) then
        print *, 'delta =', abs(z-z1)
        stop 'too many iterations in gauss_hermite';
    end if
    x(i) = z
    x(n+1-i) = -z
    w(i) = real(2.0, kind=eprec)/(pp*pp)
    w(n+1-i) = w(i)
    end do


    deallocate( d )
    deallocate( e )
    
    ! check that all nodes are distinct and increasing
    do i=2,n
        if(.not.((x(i)-x(i-1))>EPS1)) then
            stop 'nodes not distinct in gauss_hermite'
        end if
    end do    
end subroutine gauss_hermite


subroutine gauss_hermite_scaled(x, w, n, gg)
    real(kind=eprec), intent(out) :: x(n)
    real(kind=eprec), intent(out) :: w(n)
    integer, intent(in) :: n
    real(kind=eprec), intent(in) :: gg

    real(kind=eprec) :: one_over_sqrt_gg 
    integer :: m, i

    one_over_sqrt_gg = 1.0_eprec/sqrt(gg)

    call gauss_hermite(x, w, n)

    m = (n+1)/2
    do i=1,m
        w(i) = exp(x(i)**2) * w(i) * one_over_sqrt_gg
        x(i) = x(i) * one_over_sqrt_gg
        w(n+1-i) = w(i)
        x(n+1-i) = -x(i)
    end do

end subroutine gauss_hermite_scaled

subroutine hermite_polynomial(x, n, H)
    real(kind=eprec), intent(in) :: x 
    integer, intent(in) :: n
    real(kind=eprec), intent(out) :: H(0:n)

    real(kind=eprec) :: two_times_x
    integer :: i

    if (n>=0) then
       H(0) = 1.0_eprec
    end if   
    if (n>=1) then
       two_times_x = 2.0_eprec*x
       H(1) = two_times_x
       do i=2,n
           H(i) = two_times_x*H(i-1) - real(2*(i-1),kind=eprec)*H(i-2)
       end do
    end if   
end subroutine hermite_polynomial

#ifdef _QUADPRECISION_
subroutine hermite_scaled_coeffs(n, x, w, H, gg)
#else
subroutine hermite_scaled_coeffs_eprec(n, x, w, H, gg)
#endif
    integer, intent(in) :: n
    real(kind=eprec), intent(out) :: x(0:n)
    real(kind=eprec), intent(out) :: w(0:n)
    real(kind=eprec), intent(out) :: H(0:n, 0:n)
    real(kind=eprec), intent(in) :: gg 

    real(kind=eprec), parameter :: pi = 3.1415926535897932384626433832795028841971693993751_eprec
    real(kind=eprec) :: sqrt_gg 
    real(kind=eprec) :: sqrt4_gg_over_pi 
    integer :: l,k
    real(kind=eprec) :: hh(0:n)
    real(kind=eprec) :: f_l(0:n)
    real(kind=eprec) :: f_k 
    real(kind=eprec) :: u 

    sqrt_gg = sqrt(gg)
    sqrt4_gg_over_pi = sqrt(sqrt(gg/pi))

    call gauss_hermite_scaled(x, w, n+1, gg)

    u = 1.0_eprec
    f_l(0) = 1.0_eprec
    do l=1,n
        u = 2.0_eprec*u*real(l, kind=eprec)
        f_l(l) = 1.0_eprec/sqrt(u)
    end do

    do k=0,n
        f_k = exp(-0.5_eprec*gg*x(k)**2)*sqrt4_gg_over_pi
        call hermite_polynomial(sqrt_gg*x(k), n, hh)
        do l=0,n
            H(k,l) = f_k*f_l(l)*hh(l)
        end do
    end do    

#ifdef _QUADPRECISION_
end subroutine hermite_scaled_coeffs
#else
end subroutine hermite_scaled_coeffs_eprec
#endif



#ifndef _QUADPRECISION_
subroutine hermite_scaled_coeffs(n, x, w, H, gg)
    integer, intent(in) :: n
    real(kind=prec), intent(out) :: x(0:n)
    real(kind=prec), intent(out) :: w(0:n)
    real(kind=prec), intent(out) :: H(0:n, 0:n)
    real(kind=prec), intent(in) :: gg 

    real(eprec) :: x_eprec(0:n)
    real(eprec) :: w_eprec(0:n)
    real(eprec), allocatable :: H_eprec(:,:)

    allocate (H_eprec(0:n, 0:n))

    call hermite_scaled_coeffs_eprec(n, x_eprec, w_eprec, H_eprec, real(gg, kind=eprec))

    x = real(x_eprec, kind=prec) ! conversion eprec->double should be done implicitely!!!
    w = real(w_eprec, kind=prec)
    H = real(H_eprec, kind=prec)

    deallocate (H_eprec)

end subroutine hermite_scaled_coeffs
#endif



end module tssm_hermite_common
