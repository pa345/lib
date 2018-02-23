lib_LTLIBRARIES = libtrack.la

libtrack_la_CFLAGS = -fopenmp
libtrack_la_SOURCES = track.c track_filter.c track_flag.c track_offsets.c track_synth.c track_weight.c

check_PROGRAMS = print print_sc stage1

print_SOURCES = print.c
print_CFLAGS = -fopenmp
print_LDFLAGS = -fopenmp
print_LDADD = libtrack.la $(top_builddir)/curvefit/libcurvefit.la $(top_builddir)/lls/liblls.la $(top_builddir)/pomme/libpomme.la $(top_builddir)/estist/libestist_calc.la -lapex -lmagcoord -lcommon -lmsynth -lm -lcdf -lsatdata -lindices ~/usr/lib/libgsl.a -lgslcblas -lgfortran

print_sc_SOURCES = print_sc.c
print_sc_LDADD = libtrack.la -lcommon -lmsynth -lm -lcdf -lsatdata -lindices ~/usr/lib/libgsl.a -lgslcblas

stage1_SOURCES = stage1.c
stage1_LDFLAGS = -fopenmp
stage1_LDADD = libtrack.la $(top_builddir)/curvefit/libcurvefit.la $(top_builddir)/pomme/libpomme.la $(top_builddir)/estist/libestist_calc.la -lapex -lmagcoord -lcommon -lmsynth -lm -lcdf -lmwmclmcrrt -ldifi -lmldifi -lchaos -lmlchaos -lsatdata -lindices -lgfortran ~/usr/lib/libgsl.a -lgslcblas

AM_CPPFLAGS = -I$(top_builddir)/curvefit -I$(top_builddir)/chaos -I$(top_builddir)/pomme -I$(top_builddir)/estist
