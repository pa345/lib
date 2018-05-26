lib_LTLIBRARIES = libpej.la

libpej_la_SOURCES = pej.c

bin_PROGRAMS = stage1 stage2

commonlibs = $(top_builddir)/magdata/libmagdata.la $(top_builddir)/track/libtrack.la $(top_builddir)/curvefit/libcurvefit.la $(top_builddir)/lls/liblls.la $(top_builddir)/magfit/libmagfit.la $(top_builddir)/green/libgreen.la $(top_builddir)/pca/libpca.la  $(top_builddir)/estist/libestist_calc.la -lobsdata -lapex -lpomme -lmagcoord -lcommon -lmsynth -lm -lcdf -lsatdata -lindices ~/usr/lib/libgsl.a -lgslcblas -lgfortran

stage1_SOURCES = stage1.c
stage1_LDADD = ${commonlibs}

stage2_SOURCES = stage2.c
stage2_LDADD = ${commonlibs}

AM_CPPFLAGS = -I$(top_builddir)/curvefit -I$(top_builddir)/chaos -I$(top_builddir)/pomme -I$(top_builddir)/estist -I$(top_builddir)/track -I$(top_builddir)/magdata -I$(top_builddir)/magfit
