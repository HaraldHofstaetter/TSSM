#ifdef _QUADPRECISION_

module tssmq_polar
    use tssmq_polar_2D
    use tssmq_cylindrical_3D
    use tssmq_polar_real_2D
    use tssmq_cylindrical_real_3D
end module tssmq_polar

#else

module tssm_polar
    use tssm_polar_2D
    use tssm_cylindrical_3D
    use tssm_polar_real_2D
    use tssm_cylindrical_real_3D
end module tssm_polar

#endif
