lib_LTLIBRARIES = libpca3d.la

libpca3d_la_SOURCES = io.c tiegcm3d_alloc.c tiegcm3d_print.c tiegcm3d_read.c window.c

AM_CPPFLAGS = -I$(top_builddir)/green -I$(top_builddir)/magdata -I$(top_builddir)/poltor -I$(top_builddir)/lls -I$(top_builddir)/lapack_wrapper -I$(top_builddir)/magfit

common_libs = -L/usr/lib64 -lsatdata -lindices ~/usr/lib/libgsl.a /home/palken/usr/lib/liblapacke.a /home/palken/usr/lib/liblapack.a -lptf77blas -lptcblas -latlas -lpthread /usr/lib64/libfftw3.so.3 -lgfortran -lnetcdf -lcdf

bin_PROGRAMS = print_modes print_modes3b print print_fft print_fft2b spectrogram spectrogram2b stage2a stage3a stage1b stage2b stage3b print_SH print_magfield

print_SOURCES = print.c
print_LDADD = libpca3d.la -lcommon -lm -lgsl -lnetcdf ${common_libs}

print_magfield_SOURCES = print_magfield.c
print_magfield_CFLAGS = -fopenmp
print_magfield_LDFLAGS = -fopenmp
print_magfield_LDADD = libpca3d.la $(top_builddir)/magfit/libmagfit.la $(top_builddir)/pca/libpca.la $(top_builddir)/green/libgreen.la  $(top_builddir)/mageq/libmageq.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la  -ltrack -lcommon -lm -lmagfield -lmsynth -lbspline2 -lgsl -lnetcdf ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

print_fft_SOURCES = print_fft.c
print_fft_CFLAGS = -fopenmp
print_fft_LDFLAGS = -fopenmp
print_fft_LDADD = libpca3d.la -lcommon -lm -lmagfield -lgsl -lshtns_omp -lfftw3 -lfftw3_omp -lnetcdf ${common_libs}

print_fft2b_SOURCES = print_fft2b.c
print_fft2b_CFLAGS = -fopenmp
print_fft2b_LDFLAGS = -fopenmp
print_fft2b_LDADD = libpca3d.la -lcommon -lm -lmagfield -lgsl -lshtns_omp -lfftw3 -lfftw3_omp -lnetcdf ${common_libs}

print_modes_SOURCES = print_modes.c
print_modes_CFLAGS = -fopenmp
print_modes_LDFLAGS = -fopenmp
print_modes_LDADD = libpca3d.la -lcommon -lm -lmagfield -lgsl -lshtns_omp -lfftw3 -lfftw3_omp -lnetcdf ${common_libs}

print_modes3b_SOURCES = print_modes3b.c
print_modes3b_CFLAGS = -fopenmp
print_modes3b_LDFLAGS = -fopenmp
print_modes3b_LDADD = libpca3d.la -lcommon -lm -lmagfield -lgsl -lshtns_omp -lfftw3 -lfftw3_omp -lnetcdf ${common_libs}

spectrogram_SOURCES = spectrogram.c
spectrogram_LDADD = libpca3d.la -lcommon -lm -lgsl -lnetcdf ${common_libs}

spectrogram2b_SOURCES = spectrogram2b.c
spectrogram2b_CFLAGS = -fopenmp
spectrogram2b_LDFLAGS = -fopenmp
spectrogram2b_LDADD = libpca3d.la -lcommon -lm -lmagfield -lgsl -lshtns_omp -lfftw3 -lfftw3_omp -lnetcdf ${common_libs}

stage2a_SOURCES = stage2a.c
stage2a_CFLAGS = -fopenmp
stage2a_LDFLAGS = -fopenmp
stage2a_LDADD = libpca3d.la -lcommon -lm -lmagfield -lgsl -lshtns_omp -lfftw3 -lfftw3_omp -lnetcdf ${common_libs}

stage3a_SOURCES = stage3a.c
stage3a_CFLAGS = -fopenmp
stage3a_LDFLAGS = -fopenmp
stage3a_LDADD = libpca3d.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lcommon -lm -lmagfield -lgsl -lshtns_omp -lfftw3 -lfftw3_omp -lnetcdf ${common_libs}

stage1b_SOURCES = stage1b.c
stage1b_CFLAGS = -fopenmp
stage1b_LDFLAGS = -fopenmp
stage1b_LDADD = libpca3d.la -lcommon -lm -lmagfield ~/usr/lib/libgsl.a -lptcblas -latlas -lpthread -lshtns_omp -lfftw3 -lfftw3_omp -lnetcdf

stage2b_SOURCES = stage2b.c
stage2b_CFLAGS = -fopenmp
stage2b_LDFLAGS = -fopenmp
stage2b_LDADD = libpca3d.la -lcommon -lm -lmagfield ~/usr/lib/libgsl.a -lptcblas -latlas -lpthread -lshtns_omp -lfftw3 -lfftw3_omp -lnetcdf

stage3b_SOURCES = stage3b.c
stage3b_CFLAGS = -fopenmp
stage3b_LDFLAGS = -fopenmp
stage3b_LDADD = libpca3d.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la -lcommon -lm -lmagfield -lgsl -lshtns_omp -lfftw3 -lfftw3_omp -lnetcdf ${common_libs}

print_SH_SOURCES = print_SH.c
print_SH_CFLAGS = -fopenmp
print_SH_LDFLAGS = -fopenmp
print_SH_LDADD = libpca3d.la -lcommon -lm -lmagfield ~/usr/lib/libgsl.a -lptcblas -latlas -lpthread -lshtns_omp -lfftw3 -lfftw3_omp -lnetcdf

