lib_LTLIBRARIES = libmagfit.la

libmagfit_la_SOURCES = gauss.c gaussint.c magfit.c secs1d.c secs2d.c pca.c rc.c sheet.c

AM_CPPFLAGS = -I$(top_builddir)/green -I$(top_builddir)/lapack_wrapper -I$(top_builddir)/mageq -I$(top_builddir)/grobs -I$(top_builddir)/pca

check_PROGRAMS = compare main

compare_SOURCES = compare.c
compare_LDADD = libmagfit.la $(top_builddir)/grobs/libgrobs.la -lcommon -lm -lindices ~/usr/lib/libgsl.a -lgslcblas

main_SOURCES = main.c
main_LDADD = libmagfit.la $(top_builddir)/pca/libpca.la $(top_builddir)/green/libgreen.la  $(top_builddir)/mageq/libmageq.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -ltrack -lcommon -lmsynth -lm -lcdf -lsatdata -lindices ~/usr/lib/libgsl.a -llapacke -llapack -lptcblas -lptf77blas -latlas -lpthread -lgfortran
