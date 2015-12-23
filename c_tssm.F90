module tssm_c_interface
    use tssm
    implicit none
contains
    subroutine c_initialize_tssm() bind(c)
        call initialize_tssm
    end subroutine c_initialize_tssm 

end module tssm_c_interface
