check_PROGRAMS = main

main_SOURCES = main.c
main_LDADD = $(top_builddir)/green/libgreen.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lm -lcommon /home/palken/usr/lib/libgsl.a /home/palken/usr/lib/liblapacke.a /home/palken/usr/lib/liblapack.a -lptf77blas -lptcblas -latlas -lgfortran -lpthread

AM_CPPFLAGS = -I$(top_builddir)/green -I$(top_builddir)/lapack_wrapper
