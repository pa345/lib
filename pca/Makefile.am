lib_LTLIBRARIES = libpca.la

libpca_la_SOURCES = io.c pca.c tiegcm_alloc.c tiegcm_read.c window.c

AM_CPPFLAGS = -I$(top_builddir)/green -I$(top_builddir)/poltor -I$(top_builddir)/lls -I$(top_builddir)/lapack_wrapper

common_libs = -L/usr/lib64 -lmainlib -lbspline2 -lsofa_c -lcspice -lmwmclmcrrt -lmlchaos ~/usr/lib/libgsl.a -lm -lcdf -lgfortran -lfftw3 -lz -lconfig /home/palken/usr/lib/liblapacke.a /home/palken/usr/lib/liblapack.a -lptf77blas -lptcblas -latlas -lpthread -lnetcdf

bin_PROGRAMS = plot_maps print print2 stage1 stage2 stage2b stage3 stage3b stage1_ind stage2b_ind

plot_maps_SOURCES = plot_maps.c
plot_maps_LDADD = libpca.la $(top_builddir)/green/libgreen.la ${common_libs}

print_SOURCES = print.c
print_LDADD = libpca.la $(top_builddir)/green/libgreen.la ${common_libs}

print2_SOURCES = print2.c
print2_LDADD = libpca.la $(top_builddir)/green/libgreen.la ${common_libs}

stage1_SOURCES = stage1.c
stage1_CFLAGS = -fopenmp
stage1_LDFLAGS = -fopenmp
stage1_LDADD = libpca.la $(top_builddir)/lls/liblls.la $(top_builddir)/green/libgreen.la  $(top_builddir)/lapack_wrapper/liblapack_wrapper.la ${common_libs}

stage1_ind_SOURCES = stage1_ind.c
stage1_ind_CFLAGS = -fopenmp
stage1_ind_LDFLAGS = -fopenmp
stage1_ind_LDADD = libpca.la $(top_builddir)/lls/liblls.la $(top_builddir)/green/libgreen.la  $(top_builddir)/lapack_wrapper/liblapack_wrapper.la ${common_libs}

stage2_SOURCES = stage2.c
stage2_LDADD = libpca.la $(top_builddir)/green/libgreen.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la ${common_libs}

stage2b_SOURCES = stage2b.c
stage2b_LDADD = libpca.la $(top_builddir)/green/libgreen.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la ${common_libs}

stage2b_ind_SOURCES = stage2b_ind.c
stage2b_ind_LDADD = libpca.la $(top_builddir)/green/libgreen.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la ${common_libs}

stage3_SOURCES = stage3.c
stage3_LDADD = libpca.la $(top_builddir)/green/libgreen.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la ${common_libs}

stage3b_SOURCES = stage3b.c
stage3b_LDADD = libpca.la $(top_builddir)/green/libgreen.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la ${common_libs}
