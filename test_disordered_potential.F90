program test_disorder_potential
    use tssm_schroedinger
    use tssm_disorder_potential
    implicit none

    type(schroedinger_1D) :: m
    real(prec) :: mean, dev, V0, mu, xi

    call initialize_tssm


    mu = 0.349822723719617E+02_prec
    xi = 1.0_prec/sqrt(2.0_prec*mu)

    m = schroedinger_1D(nx = 2048,             &
                       xmin = -400.0_prec,       &
                       xmax =  400.0_prec,       &
!                       ny = 2048,              &
!                       ymin = -400.0_prec,       &
!                       ymax =  400.0_prec,       &
                       hbar = 1.0_prec,          &
                       mass = 1.0_prec,          &
                       cubic_coupling = 390.0_prec, &
!                       cubic_coupling = 5274.0_prec, &
                       potential = weak_periodic_1D,    &
                       boundary_conditions = periodic)

    call compute_disorder_potential_1D(m%g, m%V, &
                   dx = 1.0_prec,      &
!                   dy = 0.1_prec,      &
                   sigma = 0.7_prec*xi, &
                   f = 10)        


   call get_mean_value_and_standard_deviation_1D(m%g, m%V, mean, dev)
   print *, "mean, dev =", mean, dev

   call m%save_potential("V_disorder_1D.h5")
   !call m%load_potential("V_disorder_1D.h5")

   call get_mean_value_and_standard_deviation_1D(m%g, m%V, mean, dev)
   print *, "mean, dev =", mean, dev

   !normalize potential
   V0 = 1.0_prec
   m%V = (V0/dev)*(m%V - mean)

   call get_mean_value_and_standard_deviation_1D(m%g, m%V, mean, dev)
   print *, "mean, dev =", mean, dev

   !call m%save_potential("V_disorder_1D.h5")


    call finalize_tssm

contains

    function weak_periodic_1D(x) result(y)
        real(kind=prec), intent(in) :: x
        real(kind=prec) :: y
        y = 1.4_prec * cos(11.46_prec*x) 
    end function weak_periodic_1D

    function weak_periodic_2D(x,y) result(z)
        real(kind=prec), intent(in) :: x,y
        real(kind=prec) :: z
        z = 1.4_prec * cos(11.46_prec*x) * cos(11.46*y)
    end function weak_periodic_2D


end program test_disorder_potential
