#ifdef _QUADPRECISION_

module tssmq_base
     use tssmq_common
     use tssmq_base_complex
     use tssmq_base_real
end module tssmq_base

#else

module tssm_base
     use tssm_common
     use tssm_base_complex
     use tssm_base_real
end module tssm_base

#endif


