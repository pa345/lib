lib_LTLIBRARIES = libephemeris.la
libephemeris_la_SOURCES = eph.c eph_data.c hermite.c

AM_CPPFLAGS =

bin_PROGRAMS = print

print_SOURCES = print.c
print_LDADD = libephemeris.la -lcommon -lm -lcdf -lsatdata -lindices -lz -lgsl -lgslcblas

check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libephemeris.la -lcommon -lm -lcdf -lsatdata -lindices -lz -lsofa -lsofa_c -lspice -lcspice -lgsl -lgslcblas
