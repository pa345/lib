check_PROGRAMS = necorr necorrJ

necorr_SOURCES = necorr.c
necorr_LDADD = $(top_builddir)/peak/libpeak.la $(top_builddir)/euler/libeuler.la -ltrack -lcommon -lmsynth -lm -lcdf -lsatdata -lindices ~/usr/lib/libgsl.a -lgslcblas

necorrJ_SOURCES = necorrJ.c
necorrJ_LDADD = $(top_builddir)/peak/libpeak.la $(top_builddir)/euler/libeuler.la -ltrack -lapex -lcommon -lmsynth -lm -lcdf -lsatdata -lindices ~/usr/lib/libgsl.a -lgslcblas -lgfortran

AM_CPPFLAGS = -I$(top_builddir)/peak
