#ifdef _QUADPRECISION_

module tssmq_tensorial
    use tssmq_tensorial_1D
    use tssmq_tensorial_2D
    use tssmq_tensorial_3D
    use tssmq_tensorial_real_1D
    use tssmq_tensorial_real_2D
    use tssmq_tensorial_real_3D
end module tssmq_tensorial

#else

module tssm_tensorial
    use tssm_tensorial_1D
    use tssm_tensorial_2D
    use tssm_tensorial_3D
    use tssm_tensorial_real_1D
    use tssm_tensorial_real_2D
    use tssm_tensorial_real_3D
end module tssm_tensorial

#endif
