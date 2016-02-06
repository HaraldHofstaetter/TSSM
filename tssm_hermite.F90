#ifdef _QUADPRECISION_

module tssmq_hermite
    use tssmq_hermite_1D
    use tssmq_hermite_2D
    use tssmq_hermite_3D
    use tssmq_hermite_real_1D
    use tssmq_hermite_real_2D
    use tssmq_hermite_real_3D
end module tssmq_hermite

#else

module tssm_hermite
    use tssm_hermite_1D
    use tssm_hermite_2D
    use tssm_hermite_3D
    use tssm_hermite_real_1D
    use tssm_hermite_real_2D
    use tssm_hermite_real_3D
end module tssm_hermite

#endif

