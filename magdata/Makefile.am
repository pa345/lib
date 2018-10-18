noinst_LTLIBRARIES = libmagdata.la

AM_CPPFLAGS =

libmagdata_la_SOURCES = magdata.c magdata_list.c preproc.c

check_PROGRAMS = test print

print_SOURCES = print.c
print_LDADD = libmagdata.la -ltrack -lobsdata -lapex -lcommon -lmsynth -lm -lgsl -lgslcblas -lsatdata -lindices -lpomme -lcdf -lgfortran -L/home/palken/usr/lib

test_SOURCES = test.c
test_LDADD = libmagdata.la -ltrack -lobsdata -lapex -lcommon -lmsynth -lm -lgsl -lgslcblas -lsatdata -lindices -lpomme -lcdf -lgfortran -L/home/palken/usr/lib
