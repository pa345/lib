lib_LTLIBRARIES = libtrack.la

libtrack_la_CFLAGS = -fopenmp
libtrack_la_SOURCES = track.c track_flag.c track_synth.c track_weight.c

check_PROGRAMS = print print_sc stage1

print_SOURCES = print.c
print_CFLAGS = -fopenmp
print_LDFLAGS = -fopenmp
print_LDADD = libtrack.la -lpomme -lapex -lmagcoord -lcommon -lmsynth -lbspline2 -lm -lcdf -lsatdata -lindices ~/usr/lib/libgsl.a -lgslcblas -lgfortran

print_sc_SOURCES = print_sc.c
print_sc_LDADD = libtrack.la -lcommon -lmsynth -lbspline2 -lm -lcdf -lsatdata -lindices ~/usr/lib/libgsl.a -lgslcblas

stage1_SOURCES = stage1.c
stage1_LDFLAGS = -fopenmp
stage1_LDADD = libtrack.la -lpomme -lapex -lmagcoord -lcommon -lmsynth -lbspline2 -lm -lcdf -lmwmclmcrrt -ldifi -lmldifi -lchaos -lmlchaos -lsatdata -lindices -lgfortran ~/usr/lib/libgsl.a -lgslcblas

AM_CPPFLAGS =
