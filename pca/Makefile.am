lib_LTLIBRARIES = libpca.la

libpca_la_SOURCES = io.c pca.c tiegcm_alloc.c tiegcm_read.c window.c

AM_CPPFLAGS = -I$(top_builddir)/green -I$(top_builddir)/magdata -I$(top_builddir)/track -I$(top_builddir)/poltor -I$(top_builddir)/lls -I$(top_builddir)/lapack_wrapper

common_libs = $(top_builddir)/curvefit/libcurvefit.la $(top_builddir)/pomme/libpomme.la $(top_builddir)/estist/libestist_calc.la -L/usr/lib64 -lobsdata -lmagcoord -lapex -lpomme -lcommon -lmsynth -lsatdata -lindices ~/usr/lib/libgsl.a /home/palken/usr/lib/liblapacke.a /home/palken/usr/lib/liblapack.a -lptf77blas -lptcblas -latlas -lpthread /usr/lib64/libfftw3.so.3 -lgfortran -lnetcdf -lcdf

bin_PROGRAMS = plot_maps print print2 stage1 stage2 stage2b stage3 stage3b stage1_ind stage2b_ind

plot_maps_SOURCES = plot_maps.c
plot_maps_LDADD = libpca.la $(top_builddir)/green/libgreen.la -lcommon -lm ${common_libs}

print_SOURCES = print.c
print_LDADD = libpca.la $(top_builddir)/green/libgreen.la -lcommon -lm -lnetcdf ${common_libs}

print2_SOURCES = print2.c
print2_LDADD = libpca.la $(top_builddir)/green/libgreen.la -lcommon -lm -lnetcdf ${common_libs}

stage1_SOURCES = stage1.c
stage1_CFLAGS = -fopenmp
stage1_LDFLAGS = -fopenmp
stage1_LDADD = libpca.la $(top_builddir)/magdata/libmagdata.la $(top_builddir)/lls/liblls.la $(top_builddir)/green/libgreen.la  $(top_builddir)/track/libtrack.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lm -lcdf -lnetcdf ${common_libs}

stage1_ind_SOURCES = stage1_ind.c
stage1_ind_CFLAGS = -fopenmp
stage1_ind_LDFLAGS = -fopenmp
stage1_ind_LDADD = libpca.la $(top_builddir)/magdata/libmagdata.la $(top_builddir)/lls/liblls.la $(top_builddir)/green/libgreen.la  $(top_builddir)/track/libtrack.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lm -lcdf -lnetcdf ${common_libs}

stage2_SOURCES = stage2.c
stage2_LDADD = libpca.la $(top_builddir)/green/libgreen.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lcommon -lm ${common_libs}

stage2b_SOURCES = stage2b.c
stage2b_LDADD = libpca.la $(top_builddir)/green/libgreen.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lcommon -lm ${common_libs}

stage2b_ind_SOURCES = stage2b_ind.c
stage2b_ind_LDADD = libpca.la $(top_builddir)/green/libgreen.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lcommon -lm ${common_libs}

stage3_SOURCES = stage3.c
stage3_LDADD = libpca.la $(top_builddir)/green/libgreen.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lcommon -lm ${common_libs}

stage3b_SOURCES = stage3b.c
stage3b_LDADD = libpca.la $(top_builddir)/green/libgreen.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lcommon -lm ${common_libs}
