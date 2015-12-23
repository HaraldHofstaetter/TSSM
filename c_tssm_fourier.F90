module c_tssm_fourier
    use tssm_fourier
    implicit none
contains

    subroutine c_initialize_tssm_fourier() bind(c)
        call initialize_tssm_fourier
    end subroutine c_initialize_tssm_fourier 

end module c_tssm_fourier
