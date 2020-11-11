lib_LTLIBRARIES = libmagfit.la

libmagfit_la_SOURCES = gauss.c gaussint.c magfit.c lc_scalar.c secs1d.c secs2d.c pca.c rc.c sheet.c

AM_CPPFLAGS = -I$(top_builddir)/grobs -I$(top_builddir)/pca

check_PROGRAMS = compare main

common_libs = -lmainlib -lbspline2 -lsofa_c -lcspice -lmwmclmcrrt -lmlchaos ~/usr/lib/libgsl.a -lm -lcdf -lgfortran -lfftw3 -lz -lconfig -lptcblas -lptf77blas -latlas -lpthread

compare_SOURCES = compare.c
compare_LDADD = libmagfit.la $(top_builddir)/grobs/libgrobs.la ${common_libs}

main_SOURCES = main.c
main_LDADD = libmagfit.la $(top_builddir)/pca/libpca.la -llapacke -llapack ${common_libs}
