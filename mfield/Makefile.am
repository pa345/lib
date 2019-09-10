bin_PROGRAMS = mfield mfield_preproc mfield_preproc_EEJ mfield_preproc_obs

check_PROGRAMS = mfield_emag mfield_eval_main mfield_compare mfield_plot

common_libs = -L/home/palken/usr/lib -lmainlib -lbspline2 -lsofa_c -lcspice -lmwmclmcrrt -lmlchaos ~/usr/lib/libgsl.a -lm -lcdf -lgfortran -lfftw3 -lz -lconfig -lflow -lm -lgslspblas2 ~/usr/lib/liblapacke.a ~/usr/lib/liblapack.a ~/usr/lib/libptcblas.a ~/usr/lib/libptf77blas.a ~/usr/lib/libatlas.a -lpthread

mfield_SOURCES = mfield.c mfield_data.c mfield_error.c mfield_align.c mfield_fluxcal.c mfield_main.c mfield_residual.c mfield_synth.c
mfield_CFLAGS = -fopenmp
mfield_LDFLAGS = -fopenmp
mfield_LDADD = ${common_libs}

mfield_emag_SOURCES = mfield_emag.c
mfield_emag_LDADD = ${common_libs}

mfield_preproc_SOURCES = mfield_preproc.c
mfield_preproc_CFLAGS = -fopenmp
mfield_preproc_LDFLAGS = -fopenmp
mfield_preproc_LDADD = ${common_libs}

mfield_preproc_EEJ_SOURCES = mfield_preproc_EEJ.c
mfield_preproc_EEJ_CFLAGS = -fopenmp
mfield_preproc_EEJ_LDFLAGS = -fopenmp
mfield_preproc_EEJ_LDADD = ${common_libs}

mfield_preproc_obs_SOURCES = mfield_preproc_obs.c
mfield_preproc_obs_CFLAGS = -fopenmp
mfield_preproc_obs_LDFLAGS = -fopenmp
mfield_preproc_obs_LDADD = ${common_libs}

mfield_compare_SOURCES = mfield_compare.c
mfield_compare_LDADD = ${common_libs}

mfield_plot_SOURCES = mfield_plot.c
mfield_plot_LDADD = ${common_libs}

mfield_eval_main_SOURCES = mfield_eval_main.c mfield_eval.c
mfield_eval_main_LDADD = -lm /home/palken/usr/lib/libgsl.a -lgslcblas
