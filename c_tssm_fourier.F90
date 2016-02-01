module c_tssm_fourier
    use tssm_fourier
    implicit none
contains

    subroutine c_initialize_tssm_fourier() bind(c)
        call initialize_tssm_fourier
    end subroutine c_initialize_tssm_fourier 

    subroutine c_set_fftw_planning_rigor(flag) bind(c)
        integer(c_int), value :: flag
        fftw_planning_rigor = flag
    end subroutine c_set_fftw_planning_rigor

end module c_tssm_fourier
