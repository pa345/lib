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

AM_CPPFLAGS = -I/home/palken/usr/include -I$(top_builddir)/cond -I$(top_builddir)/hwm -I$(top_builddir)/iri -I$(top_builddir)/msis -I$(top_builddir)/lis -I$(top_builddir)/julia -I$(top_builddir)/efi -I$(top_builddir)/superlu -I$(top_builddir)/estist -I$(top_builddir)/magfit

common_libs = -lmainlib -lbspline2 -lsofa_c -lcspice -lmwmclmcrrt -lmlchaos ~/usr/lib/libgsl.a -lm -lcdf -lgfortran -lfftw3 -lz -lconfig -llapacke -llapack -lptcblas -lptf77blas -latlas -lpthread

main_LDADD = libmag.la $(top_builddir)/magfit/libmagfit.la $(top_builddir)/pca/libpca.la $(top_builddir)/hwm/libhwm.la $(top_builddir)/cond/libcond.la $(top_builddir)/iri/libiri.la $(top_builddir)/msis/libmsis.la $(top_builddir)/superlu/libsuperlu.la $(top_builddir)/estist/libestist_calc.la -lgslpde -lsuperlu_mt_OPENMP -lmagfield ${common_libs} -lshtns -lfftw3 -lfftw3_omp -lnetcdf

swarmeef_print_LDADD = libmag.la ${common_libs}

compare_LDADD = $(top_builddir)/julia/libjulia.la $(top_builddir)/efi/libefi.la ${common_libs}

corr_LDADD = ${common_libs}
