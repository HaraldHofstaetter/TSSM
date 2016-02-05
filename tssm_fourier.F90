#ifdef _QUADPRECISION_

module tssmq_fourier
    use tssmq_fourier_1D
    use tssmq_fourier_2D
    use tssmq_fourier_3D
    use tssmq_fourier_real_1D
    use tssmq_fourier_real_2D
    use tssmq_fourier_real_3D
    use tssmq_multicomponent_fourier_1D
    use tssmq_multicomponent_fourier_2D
    use tssmq_multicomponent_fourier_3D
    use tssmq_multicomponent_fourier_real_1D
    use tssmq_multicomponent_fourier_real_2D
    use tssmq_multicomponent_fourier_real_3D
end module tssmq_fourier

#else

module tssm_fourier
    use tssm_fourier_1D
    use tssm_fourier_2D
    use tssm_fourier_3D
    use tssm_fourier_real_1D
    use tssm_fourier_real_2D
    use tssm_fourier_real_3D
    use tssm_multicomponent_fourier_1D
    use tssm_multicomponent_fourier_2D
    use tssm_multicomponent_fourier_3D
    use tssm_multicomponent_fourier_real_1D
    use tssm_multicomponent_fourier_real_2D
    use tssm_multicomponent_fourier_real_3D
end module tssm_fourier

#endif


