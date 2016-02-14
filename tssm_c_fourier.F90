#ifdef _QUADPRECISION_
module tssmq_c_fourier
    use tssmq
    use tssmq_fourier
#else
module tssm_c_fourier
    use tssm
    use tssm_fourier
#endif    
    implicit none
contains

    subroutine c_initialize_tssm_fourier() &
#ifdef _QUADPRECISION_
        bind(c, name="tssmq_fourier_initialize")
#else
        bind(c, name="tssm_fourier_initialize")
#endif        
        call initialize_tssm_fourier
    end subroutine c_initialize_tssm_fourier 

    subroutine c_set_fftw_planning_rigor(flag) &
#ifdef _QUADPRECISION_
        bind(c, name="tssmq_set_fftw_planning_rigor")
#else
        bind(c, name="tssm_set_fftw_planning_rigor")
#endif        
        integer(c_int), value :: flag
        fftw_planning_rigor = flag
    end subroutine c_set_fftw_planning_rigor

#ifdef _QUADPRECISION_
end module tssmq_c_fourier
#else
end module tssm_c_fourier
#endif        
