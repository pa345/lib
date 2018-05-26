lib_LTLIBRARIES = libpoltor.la

libpoltor_la_SOURCES = poltor.c poltor_nonlinear.c poltor_shell.c poltor_synth.c

bin_PROGRAMS = invert plot preproc poltor_preproc

commonlibs = $(top_builddir)/curvefit/libcurvefit.la $(top_builddir)/track/libtrack.la $(top_builddir)/pomme/libpomme.la $(top_builddir)/estist/libestist_calc.la $(top_builddir)/green/libgreen.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -L/home/palken/usr/lib -lflow -lobsdata -lpomme -lapex -lcommon -lmsynth -lm -lcdf -lsatdata -lindices ~/usr/lib/libgsl.a -llapacke -llapack -lptcblas -lptf77blas -latlas -lpthread -lgfortran

invert_SOURCES = invert_main.c
invert_LDADD = libpoltor.la $(top_builddir)/magdata/libmagdata.la $(top_builddir)/lls/liblls.la -lm -lconfig ${commonlibs}

plot_SOURCES = plot.c
plot_LDADD = libpoltor.la $(top_builddir)/magdata/libmagdata.la $(top_builddir)/track/libtrack.la $(top_builddir)/lls/liblls.la -lm -lgfortran ${commonlibs}

preproc_SOURCES = preproc.c
preproc_LDADD = libpoltor.la $(top_builddir)/magdata/libmagdata.la $(top_builddir)/track/libtrack.la $(top_builddir)/curvefit/libcurvefit.la $(top_builddir)/lls/liblls.la $(top_builddir)/euler/libeuler.la -lm ${commonlibs}

poltor_preproc_SOURCES = poltor_preproc.c
poltor_preproc_LDADD = libpoltor.la $(top_builddir)/magdata/libmagdata.la $(top_builddir)/track/libtrack.la $(top_builddir)/curvefit/libcurvefit.la $(top_builddir)/lls/liblls.la $(top_builddir)/euler/libeuler.la -lm ${commonlibs} -lconfig

AM_CPPFLAGS = -I$(top_builddir)/pomme -I$(top_builddir)/curvefit -I$(top_builddir)/lls -I$(top_builddir)/track -I$(top_builddir)/euler -I$(top_builddir)/magdata -I$(top_builddir)/green -I$(top_builddir)/lapack_wrapper
