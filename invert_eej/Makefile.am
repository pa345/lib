lib_LTLIBRARIES = libmag.la

libmag_la_SOURCES =        \
  cfg.c                    \
	coord.c                  \
	geod2geoc.c              \
	invert.c                 \
	inverteef.c              \
  log.c                    \
	mag.c                    \
	mag_grad.c               \
	mag_eej.c                \
	mag_log.c                \
	mag_log_grad.c           \
	mag_sqfilt_scalar.c      \
	mag_sqfilt_vector.c      \
	pde.c                    \
	sigma.c                  \
	swarmeef.c

check_PROGRAMS = corr main swarmeef_print compare

main_SOURCES = main.c
swarmeef_print_SOURCES = swarmeef_print.c
compare_SOURCES = compare.c
corr_SOURCES = corr.c

AM_CFLAGS = -Wall -W

AM_CPPFLAGS = -I/home/palken/usr/include -I$(top_builddir)/mageq -I$(top_builddir)/cond -I$(top_builddir)/hwm -I$(top_builddir)/iri -I$(top_builddir)/msis -I$(top_builddir)/lis -I$(top_builddir)/julia -I$(top_builddir)/efi -I$(top_builddir)/superlu -I$(top_builddir)/estist -I$(top_builddir)/magfit -I$(top_builddir)/green -I$(top_builddir)/pomme

main_LDADD = libmag.la $(top_builddir)/magfit/libmagfit.la $(top_builddir)/pca/libpca.la $(top_builddir)/pomme/libpomme.la $(top_builddir)/green/libgreen.la $(top_builddir)/hwm/libhwm.la $(top_builddir)/mageq/libmageq.la $(top_builddir)/cond/libcond.la $(top_builddir)/iri/libiri.la $(top_builddir)/msis/libmsis.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la $(top_builddir)/superlu/libsuperlu.la $(top_builddir)/estist/libestist_calc.la -lapex -ltrack -lcommon -lmsynth -lbspline2 -lm -lgslpde -lsatdata -lcdf -llapacke -llapack -lsuperlu_mt_OPENMP -lindices -lmagfield -lgsl -lptcblas -lptf77blas -latlas -lgfortran -lpthread -lconfig -lshtns_omp -lfftw3 -lfftw3_omp

swarmeef_print_LDADD = libmag.la -lcommon -lsatdata -lcdf -lindices -lgsl -lgslcblas -lm

compare_LDADD = $(top_builddir)/julia/libjulia.la $(top_builddir)/efi/libefi.la -lcommon -lmsynth -lbspline2 -lm -lsatdata -lcdf -lindices -lcurvefit /home/palken/usr/lib/libgsl.a -lgslcblas

corr_LDADD = -lcommon -lmsynth -lbspline2 -lm -lsatdata -lcdf -lindices /home/palken/usr/lib/libgsl.a -lgslcblas
