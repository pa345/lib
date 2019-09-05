lib_LTLIBRARIES = libpej.la

libpej_la_SOURCES = pej.c

bin_PROGRAMS = stage1 stage1_obs stage2 print_obs

commonlibs = -lmainlib -lbspline2 -lsofa_c -lcspice -lmwmclmcrrt -lmlchaos -lm -lcdf -lgfortran -lfftw3 -lz -lconfig $(top_builddir)/magfit/libmagfit.la $(top_builddir)/pca/libpca.la -lnetcdf ~/usr/lib/libgsl.a /home/palken/usr/lib/liblapacke.a /home/palken/usr/lib/liblapack.a -lptf77blas -lptcblas -latlas -lpthread

stage1_SOURCES = stage1.c
stage1_CFLAGS = -fopenmp
stage1_LDFLAGS = -fopenmp
stage1_LDADD = ${commonlibs}

stage1_obs_SOURCES = stage1_obs.c
stage1_obs_CFLAGS = -fopenmp
stage1_obs_LDFLAGS = -fopenmp
stage1_obs_LDADD = ${commonlibs}

stage2_SOURCES = stage2.c
stage2_LDADD = ${commonlibs}

print_obs_SOURCES = print_obs.c
print_obs_LDADD = ${commonlibs}

AM_CPPFLAGS = -I$(top_builddir)/magfit
