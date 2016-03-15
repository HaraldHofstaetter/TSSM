#ifdef _QUADPRECISION_

module tssmq_schroedinger ! (Nonlinear) Schroedinger equation
    use tssmq_fourier_common
    use tssmq_schroedinger_1D
    use tssmq_schroedinger_2D
    use tssmq_schroedinger_3D
    use tssmq_schroedinger_real_1D
    use tssmq_schroedinger_real_2D
    use tssmq_schroedinger_real_3D
    use tssmq_schroedinger_rotating_2D
    use tssmq_schroedinger_rotating_3D
    use tssmq_schroedinger_hermite_1D
    use tssmq_schroedinger_hermite_2D
    use tssmq_schroedinger_hermite_3D
    use tssmq_schroedinger_hermite_real_1D
    use tssmq_schroedinger_hermite_real_2D
    use tssmq_schroedinger_hermite_real_3D
    use tssmq_schroedinger_gen_laguerre_2D
    use tssmq_schroedinger_gen_laguerre_hermite_3D
    use tssmq_schroedinger_gen_laguerre_real_2D
    use tssmq_schroedinger_gen_laguerre_hermite_real_3D
end module tssmq_schroedinger 
#else

module tssm_schroedinger ! (Nonlinear) Schroedinger equation
    use tssm_fourier_common
    use tssm_schroedinger_1D
    use tssm_schroedinger_2D
    use tssm_schroedinger_3D
    use tssm_schroedinger_real_1D
    use tssm_schroedinger_real_2D
    use tssm_schroedinger_real_3D
    use tssm_schroedinger_rotating_2D
    use tssm_schroedinger_rotating_3D
    use tssm_schroedinger_hermite_1D
    use tssm_schroedinger_hermite_2D
    use tssm_schroedinger_hermite_3D
    use tssm_schroedinger_hermite_real_1D
    use tssm_schroedinger_hermite_real_2D
    use tssm_schroedinger_hermite_real_3D
    use tssm_schroedinger_hermite_2D
    use tssm_schroedinger_hermite_3D
    use tssm_schroedinger_hermite_real_1D
    use tssm_schroedinger_hermite_real_2D
    use tssm_schroedinger_hermite_real_3D
    use tssm_schroedinger_gen_laguerre_2D
    use tssm_schroedinger_gen_laguerre_hermite_3D
    use tssm_schroedinger_gen_laguerre_real_2D
    use tssm_schroedinger_gen_laguerre_hermite_real_3D
end module tssm_schroedinger 

#endif
 
