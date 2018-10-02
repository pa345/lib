bin_PROGRAMS = mfield mfield_preproc mfield_preproc_obs

check_PROGRAMS = mfield_emag mfield_eval_main mfield_compare mfield_plot

common_libs = $(top_builddir)/curvefit/libcurvefit.la $(top_builddir)/track/libtrack.la $(top_builddir)/estist/libestist_calc.la -L/home/palken/usr/lib -lobsdata -lmagcoord -lapex -lflow -lcommon -lmsynth -lbspline2 -lm -lcdf -lpomme -lsatdata -lindices ~/usr/lib/libgsl.a ~/usr/lib/liblapacke.a ~/usr/lib/liblapack.a ~/usr/lib/libptcblas.a ~/usr/lib/libptf77blas.a ~/usr/lib/libatlas.a -lpthread -lgfortran -lnetcdf

mfield_SOURCES = mfield.c mfield_data.c mfield_error.c mfield_green.c mfield_main.c mfield_residual.c mfield_synth.c
mfield_CFLAGS = -fopenmp
mfield_LDFLAGS = -fopenmp
mfield_LDADD = $(top_builddir)/magdata/libmagdata.la $(top_builddir)/lls/liblls.la $(top_builddir)/euler/libeuler.la $(top_builddir)/green/libgreen.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lfftw3 -lconfig -lspatwt ${common_libs}

mfield_emag_SOURCES = mfield_emag.c
mfield_emag_LDADD = $(top_builddir)/magdata/libmagdata.la $(top_builddir)/track/libtrack.la $(top_builddir)/lls/liblls.la $(top_builddir)/euler/libeuler.la -lcommon -lmsynth -lbspline2 -lm -llapack -lcdf -lsatdata -lindices -lfftw3 -lgfortran -lpthread -L/home/palken/usr/lib -lgsl -lptcblas -lptf77blas -latlas -lz

mfield_preproc_SOURCES = mfield_preproc.c
mfield_preproc_CFLAGS = -fopenmp
mfield_preproc_LDFLAGS = -fopenmp
mfield_preproc_LDADD = $(top_builddir)/magdata/libmagdata.la $(top_builddir)/green/libgreen.la $(top_builddir)/lls/liblls.la $(top_builddir)/euler/libeuler.la -lfftw3 -lconfig ${common_libs}

mfield_preproc_obs_SOURCES = mfield_preproc_obs.c
mfield_preproc_obs_CFLAGS = -fopenmp
mfield_preproc_obs_LDFLAGS = -fopenmp
mfield_preproc_obs_LDADD = $(top_builddir)/magdata/libmagdata.la $(top_builddir)/green/libgreen.la $(top_builddir)/lls/liblls.la $(top_builddir)/euler/libeuler.la -lfftw3 -lconfig ${common_libs}

mfield_compare_SOURCES = mfield_compare.c
mfield_compare_LDADD = $(top_builddir)/pomme/libpomme.la $(top_builddir)/estist/libestist_calc.la -lapex -lcommon -lmsynth -lbspline2 -lm -lgsl -lptcblas -lptf77blas -latlas -lindices -lpthread -lgfortran -L/home/palken/usr/lib

mfield_plot_SOURCES = mfield_plot.c
mfield_plot_LDADD = $(top_builddir)/pomme/libpomme.la $(top_builddir)/estist/libestist_calc.la -lapex -lcommon -lmsynth -lbspline2 -lm -lgsl -lptcblas -lptf77blas -latlas -lindices -lpthread -lgfortran -L/home/palken/usr/lib

mfield_eval_main_SOURCES = mfield_eval_main.c mfield_eval.c
mfield_eval_main_LDADD = -lm /home/palken/usr/lib/libgsl.a -lgslcblas

AM_CPPFLAGS = -I$(top_builddir)/pomme -I$(top_builddir)/mageq -I$(top_builddir)/lls -I$(top_builddir)/track -I$(top_builddir)/magdata -I$(top_builddir)/euler -I$(top_builddir)/estist -I$(top_builddir)/green -I$(top_builddir)/magfit -I$(top_builddir)/lapack_wrapper
