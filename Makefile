
ARCHFILES = Makefile Makefile.common Makefile.system \
            OPTIMIZED/Makefile DEBUG/Makefile \
            OPENMP/Makefile MPI/Makefile MPI_OPENMP/Makefile \
            QUADPRECISION/Makefile MINGW/Makefile \
            QUADPREC_OPENMP/Makefile \
            lapack_dsterf.f tssm_splitting_schemes.F90 \
            tssm_common.F90 tssm_base_CPP.F90 tssm_base.F90 tssm_grid.F90 \
            tssm_hdf5.F90 tssm_hdf5_helper.F90 \
            tssm_fourier.F90 tssm_fourier_common.F90 tssm_fourier_CPP.F90 \
            tssm_multicomponent_fourier_CPP.F90 \
            tssm_c.F90 tssm_c_fourier.F90 tssm_c_fourier_CPP.F90\
            tssm_c_fourier_f128wrapper_CPP.c \
            tssm_c_fourier_bessel_f128wrapper_CPP.c \
            tssm_c_schroedinger_f128wrapper_CPP.c \
	    propagate_b_2components_complex_2d.inc \
	    propagate_b_2components_real_2d.inc \
            tssm_tensorial.F90 tssm_tensorial_CPP.F90 \
            tssm_hermite.F90 tssm_hermite_common.F90 tssm_hermite_CPP.F90 \
            tssm_polar.F90 tssm_polar_CPP.F90 \
            tssm_generalized_laguerre.F90 tssm_generalized_laguerre_common.F90 tssm_generalized_laguerre_CPP.F90 \
            tssm_c_generalized_laguerre_CPP.F90 \
            tssm_schroedinger.F90  tssm_schroedinger_CPP.F90 \
            tssm_c_schroedinger_CPP.F90 \
            tssm_fourier_bessel_common.F90 \
            tssm_fourier_bessel_CPP.F90 \
            tssm_fourier_bessel.F90 \
            tssm_c_fourier_bessel_CPP.F90 \
            tssm_disorder_potential.F90 test_fisher.F90 \
            test.F90 test_soliton.F90 test_gaussian.F90 test_harmonic.F90 \
            test_gray_scott.F90 test_gray_scott_complex.F90 test_real.F90 \
	    test_imag_time.F90 test_gross_pitaevskii.F90 test_hermite.F90 \
            test_adaptive.F90 test_disordered_potential.F90 \
            test_coupled_nonlin_schroe.F90 \
            simulation_1D.F90 simulation_2D.F90 test_paper_fig_8.F90  \
            test_compare_space_discretizations.F90 \
            Julia/tssm_tools.jl Julia/tssm.jl Julia/tssm_fourier.jl \
            Julia/tssm_schroedinger.jl Julia/tssm_schroedinger_hermite.jl \
            Julia/tssm_generalized_laguerre.jl \
            Julia/tssm_fourier_bessel.jl \
            Julia/groundstate.jl Julia/time_stepper.jl \
            Julia/coupled_schroedinger.jl Julia/test_coupled_schroedinger.jl \
            Julia/Groundstate_demo.ipynb

tgz: tssm.tgz

tssm.tgz: $(ARCHFILES)
	tar czvf tssm.tgz $(ARCHFILES)


