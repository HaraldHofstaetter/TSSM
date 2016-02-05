#ifdef _QUADPRECISION_

module tssmq
     use tssmq_complex
     use tssmq_real
end module tssmq

#else

module tssm
     use tssm_complex
     use tssm_real
end module tssm

#endif


