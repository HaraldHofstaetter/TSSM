#ifdef _QUADPRECISION_
module tssmq_c_interface
    use tssmq
#else
module tssm_c_interface
    use tssm
#endif    

    implicit none
contains
    subroutine c_initialize() &
#ifdef _QUADPRECISION_
        bind(c, name="tssmq_initialize")
#else
        bind(c, name="tssm_initialize")
#endif        
        call initialize_tssm
    end subroutine c_initialize

#ifdef _QUADPRECISION_
end module tssmq_c_interface
#else
end module tssm_c_interface
#endif    
