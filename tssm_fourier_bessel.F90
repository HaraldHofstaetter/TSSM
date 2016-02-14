#ifdef _QUADPRECISION_

module tssmq_fourier_bessel
    use tssmq_fourier_bessel_2d
    use tssmq_fourier_bessel_real_2d
!    use tssmq_fourier_bessel_rotsym_2d
!    use tssmq_fourier_bessel_rotsym_real_2d
end module tssmq_fourier_bessel

#else

module tssm_fourier_bessel
    use tssm_fourier_bessel_2d
    use tssm_fourier_bessel_real_2d
!    use tssm_fourier_bessel_rotsym_2d
!    use tssm_fourier_bessel_rotsym_real_2d
end module tssm_fourier_bessel

#endif
