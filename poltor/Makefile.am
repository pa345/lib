lib_LTLIBRARIES = libpoltor.la

libpoltor_la_SOURCES = poltor.c poltor_nonlinear.c poltor_shell.c poltor_synth.c

bin_PROGRAMS = invert plot preproc poltor_preproc

commonlibs = -L/home/palken/usr/lib -lflow -lmainlib -lbspline2 -lsofa_c -lcspice -lmwmclmcrrt -lmlchaos ~/usr/lib/libgsl.a -lm -lcdf -lgfortran -lfftw3 -lz -lconfig -llapacke -llapack -lptcblas -lptf77blas -latlas -lpthread

invert_SOURCES = invert_main.c
invert_LDADD = libpoltor.la ${commonlibs}

plot_SOURCES = plot.c
plot_LDADD = libpoltor.la ${commonlibs}

preproc_SOURCES = preproc.c
preproc_LDADD = libpoltor.la ${commonlibs}

poltor_preproc_SOURCES = poltor_preproc.c
poltor_preproc_LDADD = libpoltor.la ${commonlibs}
