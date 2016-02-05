#ifdef _QUADPRECISION_

module tssmq_generalized_laguerre
    use tssmq_generalized_laguerre_2D
    use tssmq_generalized_laguerre_hermite_3D
    use tssmq_generalized_laguerre_real_2D
    use tssmq_generalized_laguerre_hermite_real_3D
end module tssmq_generalized_laguerre

#else

module tssm_generalized_laguerre
    use tssm_generalized_laguerre_2D
    use tssm_generalized_laguerre_hermite_3D
    use tssm_generalized_laguerre_real_2D
    use tssm_generalized_laguerre_hermite_real_3D
end module tssm_generalized_laguerre

#endif
