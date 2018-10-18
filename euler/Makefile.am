noinst_LTLIBRARIES = libeuler.la

AM_CPPFLAGS = -I$(top_builddir)/magdata

libeuler_la_SOURCES = euler.c euler_calc.c

check_PROGRAMS = test test_quat main preproc print

main_SOURCES = main.c
main_LDADD = libeuler.la $(top_builddir)/magdata/libmagdata.la -ltrack -lcommon -lmsynth -lm -lgsl -llapack -lf77blas -lcblas -latlas -lgfortran -lsatdata -lindices -lcdf -L/home/palken/usr/lib

preproc_SOURCES = preproc.c
preproc_LDADD = libeuler.la $(top_builddir)/magdata/libmagdata.la -ltrack -lcommon -lmsynth -lm -lsatdata -lindices -lcdf -L/home/palken/usr/lib -lgsl -lgslcblas

print_SOURCES = print.c
print_LDADD = libeuler.la -lcommon -lsatdata -lindices -lcdf -lm -lgsl -lgslcblas -L/home/palken/usr/lib

test_SOURCES = test.c
test_LDADD = libeuler.la -lcommon -lm -lgsl -lgslcblas -lsatdata -lindices -lcdf -L/home/palken/usr/lib

test_quat_SOURCES = test_quat.c
test_quat_LDADD = libeuler.la $(top_builddir)/magdata/libmagdata.la -ltrack -lcommon -lmsynth -lm -lgsl -lgslcblas -lsatdata -lindices -lcdf -L/home/palken/usr/lib
