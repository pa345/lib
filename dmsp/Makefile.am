bin_PROGRAMS = print print_LTAN stage0 stage2a stage2b

commonlibs = -lsatdata -lindices ~/usr/lib/libgsl.a -lgslcblas -lm -lcdf

AM_CPPFLAGS = -I$(top_builddir)/ephemeris -I$(top_builddir)/chaos -I$(top_builddir)/track -I$(top_builddir)/estist -I$(top_builddir)/pomme -I$(top_builddir)/euler -I$(top_builddir)/magdata -I$(top_builddir)/peak

stage0_SOURCES = stage0.c
stage0_LDADD = $(top_builddir)/ephemeris/libephemeris.la -lcommon -lm -lcdf ${commonlibs} -lz

stage2a_SOURCES = attitude.c jump.c fluxcal.c magcal.c stage2a.c
stage2a_LDFLAGS = -fopenmp
stage2a_LDADD = $(top_builddir)/track/libtrack.la $(top_builddir)/ephemeris/libephemeris.la $(top_builddir)/euler/libeuler.la $(top_builddir)/peak/libpeak.la -lcommon -lmsynth -lm -lgfortran ${commonlibs} -lz

stage2b_SOURCES = stage2b.c magcal.c
stage2b_LDADD = -lcommon -lm ${commonlibs}

print_SOURCES = print.c
print_LDADD = -lcommon -lm -lcdf ${commonlibs}

print_LTAN_SOURCES = print_LTAN.c
print_LTAN_LDADD = $(top_builddir)/track/libtrack.la $(top_builddir)/ephemeris/libephemeris.la $(top_builddir)/euler/libeuler.la -lcommon -lmsynth -lm -lgfortran ${commonlibs} -lz
