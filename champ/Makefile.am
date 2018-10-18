bin_PROGRAMS = print print_alt print_lt print_plp

commonlibs = ~/usr/lib/liblapack.a ~/usr/lib/libf77blas.a ~/usr/lib/libatlas.a -lcommon -lmsynth -lbspline2 -lsatdata -lindices ~/usr/lib/libgsl.a -lgslcblas -lm -lcdf -lgfortran

print_SOURCES = print.c
print_LDADD = $(top_builddir)/euler/libeuler.la -lm ${commonlibs}

print_alt_SOURCES = print_alt.c
print_alt_LDADD = -lm ${commonlibs}

print_lt_SOURCES = print_lt.c
print_lt_LDADD = -lm ${commonlibs}

print_plp_SOURCES = print_plp.c
print_plp_LDADD = -lm ${commonlibs}

AM_CPPFLAGS = -I$(top_builddir)/euler -I$(top_builddir)/magdata
