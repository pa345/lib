lib_LTLIBRARIES = libpca3d.la

libpca3d_la_SOURCES = io.c tiegcm3d_alloc.c tiegcm3d_print.c tiegcm3d_read.c window.c

AM_CPPFLAGS = -I$(top_builddir)/poltor -I$(top_builddir)/magfit

common_libs = -L/usr/lib64 -lmainlib -lbspline2 -lsofa_c -lcspice -lmwmclmcrrt -lmlchaos ~/usr/lib/libgsl.a -lm -lcdf -lgfortran -lfftw3 -lz -lconfig /home/palken/usr/lib/liblapacke.a /home/palken/usr/lib/liblapack.a -lptf77blas -lptcblas -latlas -lpthread /usr/lib64/libfftw3.so.3 -lgfortran -lnetcdf -lcdf -lcsparse

bin_PROGRAMS = invert invert_preproc invert_write_tmode print_modes print_modes3b print print_B print_fft print_fft2b spectrogram spectrogram2b stage2a stage3a stage1b stage2b stage3b print_SH1b print_map1b print_spectrum1b print_ezie_grid print_magfield print_var3b invert_print invert_postproc_lcurve invert_postproc_print

invert_SOURCES = invert_main.c invert.c invert_data.c invert_regularize.c invert_residual.c invert_tmode.c invert_smode.c invert_synth.c
#invert_CFLAGS = -std=c11 -fopenmp -ffast-math
invert_CFLAGS = -std=c11 -fopenmp
invert_LDFLAGS = -fopenmp
invert_LDADD = libpca3d.la -lmagfield -lbspline2 -lgslspblas2 ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

invert_preproc_SOURCES = invert_preproc.c
invert_preproc_CFLAGS = -fopenmp
invert_preproc_LDFLAGS = -fopenmp
invert_preproc_LDADD = libpca3d.la -lmagfield -lbspline2 ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

invert_print_SOURCES = invert_print.c invert.c invert_data.c invert_regularize.c invert_residual.c invert_tmode.c invert_smode.c invert_synth.c
invert_print_CFLAGS = -std=c11 -fopenmp
invert_print_LDFLAGS = -fopenmp
invert_print_LDADD = libpca3d.la -lmagfield -lbspline2 -lgslspblas2 ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

invert_postproc_lcurve_SOURCES = invert_postproc_lcurve.c invert.c invert_data.c invert_regularize.c invert_residual.c invert_tmode.c invert_smode.c invert_synth.c
invert_postproc_lcurve_CFLAGS = -std=c11 -fopenmp
invert_postproc_lcurve_LDFLAGS = -fopenmp
invert_postproc_lcurve_LDADD = libpca3d.la -lmagfield -lbspline2 -lgslspblas2 ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

invert_postproc_print_SOURCES = invert_postproc_print.c invert.c invert_data.c invert_regularize.c invert_residual.c invert_tmode.c invert_smode.c invert_synth.c
invert_postproc_print_CFLAGS = -std=c11 -fopenmp
invert_postproc_print_LDFLAGS = -fopenmp
invert_postproc_print_LDADD = libpca3d.la -lmagfield -lbspline2 -lgslspblas2 ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

invert_write_tmode_SOURCES = invert_write_tmode.c invert_tmode.c
invert_write_tmode_LDADD = ${common_libs}

print_SOURCES = print.c
print_LDADD = libpca3d.la -lmagfield ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

print_B_SOURCES = print_B.c
print_B_LDADD = libpca3d.la -lmagfield ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

print_magfield_SOURCES = print_magfield.c
print_magfield_CFLAGS = -fopenmp
print_magfield_LDFLAGS = -fopenmp
print_magfield_LDADD = libpca3d.la $(top_builddir)/magfit/libmagfit.la $(top_builddir)/pca/libpca.la -lmagfield -lbspline2 ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

print_ezie_grid_SOURCES = print_ezie_grid.c
print_ezie_grid_CFLAGS = -fopenmp
print_ezie_grid_LDFLAGS = -fopenmp
print_ezie_grid_LDADD = libpca3d.la $(top_builddir)/magfit/libmagfit.la $(top_builddir)/pca/libpca.la -lmagfield -lbspline2 ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

print_fft_SOURCES = print_fft.c
print_fft_CFLAGS = -fopenmp
print_fft_LDFLAGS = -fopenmp
print_fft_LDADD = libpca3d.la -lmagfield -lshtns_omp -lfftw3 -lfftw3_omp ${common_libs}

print_fft2b_SOURCES = print_fft2b.c
print_fft2b_CFLAGS = -fopenmp
print_fft2b_LDFLAGS = -fopenmp
print_fft2b_LDADD = libpca3d.la -lmagfield -lshtns_omp -lfftw3 -lfftw3_omp ${common_libs}

print_modes_SOURCES = print_modes.c
print_modes_CFLAGS = -fopenmp
print_modes_LDFLAGS = -fopenmp
print_modes_LDADD = libpca3d.la -lmagfield -lshtns_omp -lfftw3 -lfftw3_omp ${common_libs}

print_modes3b_SOURCES = print_modes3b.c
print_modes3b_CFLAGS = -fopenmp
print_modes3b_LDFLAGS = -fopenmp
print_modes3b_LDADD = libpca3d.la -lmagfield -lshtns_omp -lfftw3 -lfftw3_omp ${common_libs}

print_var3b_SOURCES = print_var3b.c
print_var3b_CFLAGS = -fopenmp
print_var3b_LDFLAGS = -fopenmp
print_var3b_LDADD = libpca3d.la -lmagfield -lshtns_omp -lfftw3 -lfftw3_omp ${common_libs}

spectrogram_SOURCES = spectrogram.c
spectrogram_LDADD = libpca3d.la -lmagfield ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

spectrogram2b_SOURCES = spectrogram2b.c
spectrogram2b_CFLAGS = -fopenmp
spectrogram2b_LDFLAGS = -fopenmp
spectrogram2b_LDADD = libpca3d.la -lmagfield -lshtns_omp -lfftw3 -lfftw3_omp ${common_libs}

stage2a_SOURCES = stage2a.c
stage2a_CFLAGS = -fopenmp
stage2a_LDFLAGS = -fopenmp
stage2a_LDADD = libpca3d.la -lmagfield -lshtns_omp -lfftw3 -lfftw3_omp ${common_libs}

stage3a_SOURCES = stage3a.c
stage3a_CFLAGS = -fopenmp
stage3a_LDFLAGS = -fopenmp
stage3a_LDADD = libpca3d.la -lmagfield -lshtns_omp -lfftw3 -lfftw3_omp ${common_libs}

stage1b_SOURCES = stage1b.c
stage1b_CFLAGS = -fopenmp
stage1b_LDFLAGS = -fopenmp
stage1b_LDADD = libpca3d.la -lmagfield ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

stage2b_SOURCES = stage2b.c
stage2b_CFLAGS = -fopenmp
stage2b_LDFLAGS = -fopenmp
stage2b_LDADD = libpca3d.la -lmagfield ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

stage3b_SOURCES = stage3b.c
stage3b_CFLAGS = -fopenmp
stage3b_LDFLAGS = -fopenmp
stage3b_LDADD = libpca3d.la -lmagfield -lshtns_omp -lfftw3 -lfftw3_omp ${common_libs}

#stage1c_SOURCES = stage1c.c
#stage1c_CFLAGS = -fopenmp
#stage1c_LDFLAGS = -fopenmp
#stage1c_LDADD = libpca3d.la -lcommon -lm -lmagfield ~/usr/lib/libgsl.a -lptcblas -latlas -lpthread -lshtns_omp -lfftw3 -lfftw3_omp -lnetcdf

print_SH1b_SOURCES = print_SH1b.c
print_SH1b_CFLAGS = -fopenmp
print_SH1b_LDFLAGS = -fopenmp
print_SH1b_LDADD = libpca3d.la -lmagfield ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

print_map1b_SOURCES = print_map1b.c
print_map1b_CFLAGS = -fopenmp
print_map1b_LDFLAGS = -fopenmp
print_map1b_LDADD = libpca3d.la -lmagfield ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

print_spectrum1b_SOURCES = print_spectrum1b.c
print_spectrum1b_CFLAGS = -fopenmp
print_spectrum1b_LDFLAGS = -fopenmp
print_spectrum1b_LDADD = libpca3d.la -lmagfield ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp
