lib_LTLIBRARIES = libpca3d.la

libpca3d_la_SOURCES = io.c tiegcm3d_alloc.c tiegcm3d_print.c tiegcm3d_read.c window.c

AM_CPPFLAGS = -I$(top_builddir)/green -I$(top_builddir)/magdata -I$(top_builddir)/poltor -I$(top_builddir)/lls -I$(top_builddir)/lapack_wrapper -I$(top_builddir)/magfit

common_libs = -L/usr/lib64 -lsatdata -lindices ~/usr/lib/libgsl.a /home/palken/usr/lib/liblapacke.a /home/palken/usr/lib/liblapack.a -lptf77blas -lptcblas -latlas -lpthread /usr/lib64/libfftw3.so.3 -lgfortran -lnetcdf -lcdf

bin_PROGRAMS = invert invert_preproc invert_write_tmode print_modes print_modes3b print print_B print_fft print_fft2b spectrogram spectrogram2b stage2a stage3a stage1b stage2b stage3b print_SH print_magfield print_fillion print_fillion2 print_gary print_var3b

invert_SOURCES = invert_main.c invert.c invert_data.c invert_residual.c invert_tmode.c invert_smode.c
#invert_CFLAGS = -std=c11 -fopenmp -ffast-math
invert_CFLAGS = -std=c11 -fopenmp
invert_LDFLAGS = -fopenmp
#invert_LDADD = libpca3d.la ${top_builddir}/magdata/libmagdata.la ${top_builddir}/green/libgreen.la -lobsdata -ltrack -lmagcoord -lapex -lconfig -lcommon -lmagfield -lmsynth -lbspline2 -lpomme -lm -lgslspblas2 ~/usr/lib/libgsl.a -lspatwt ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp
invert_LDADD = libpca3d.la ${top_builddir}/magdata/libmagdata.la ${top_builddir}/green/libgreen.la -lobsdata -ltrack -lmagcoord -lapex -lconfig -lcommon -lmagfield -lmsynth -lbspline2 -lpomme -lm -lgslspblas2 ~/usr/lib/libgsl_prof.a -lspatwt ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

invert_preproc_SOURCES = invert_preproc.c
invert_preproc_CFLAGS = -fopenmp
invert_preproc_LDFLAGS = -fopenmp
invert_preproc_LDADD = libpca3d.la ${top_builddir}/magdata/libmagdata.la -lobsdata -ltrack -lmagcoord -lapex -lconfig -lcommon -lmsynth -lbspline2 -lpomme -lm -lgsl -lnetcdf ${common_libs}

invert_write_tmode_SOURCES = invert_write_tmode.c invert_tmode.c
invert_write_tmode_LDADD = -lcommon -lm -lgsl -lnetcdf ${common_libs}

print_SOURCES = print.c
print_LDADD = libpca3d.la -lcommon -lm -lgsl -lnetcdf ${common_libs}

print_B_SOURCES = print_B.c
print_B_LDADD = libpca3d.la -lcommon -lm -lgsl -lnetcdf ${common_libs}

print_magfield_SOURCES = print_magfield.c
print_magfield_CFLAGS = -fopenmp
print_magfield_LDFLAGS = -fopenmp
print_magfield_LDADD = libpca3d.la $(top_builddir)/magfit/libmagfit.la $(top_builddir)/pca/libpca.la $(top_builddir)/green/libgreen.la  $(top_builddir)/mageq/libmageq.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la  -ltrack -lapex -lmagcoord -lpomme -lcommon -lm -lmagfield -lmsynth -lbspline2 -lgsl -lnetcdf ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

print_gary_SOURCES = print_gary.c
print_gary_CFLAGS = -fopenmp
print_gary_LDFLAGS = -fopenmp
print_gary_LDADD = libpca3d.la $(top_builddir)/green/libgreen.la  $(top_builddir)/mageq/libmageq.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la  -ltrack -lapex -lmagcoord -lpomme -lcommon -lm -lmagfield -lmsynth -lbspline2 -lgsl -lnetcdf ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

print_fillion_SOURCES = print_fillion.c
print_fillion_CFLAGS = -fopenmp
print_fillion_LDFLAGS = -fopenmp
print_fillion_LDADD = libpca3d.la $(top_builddir)/green/libgreen.la  $(top_builddir)/mageq/libmageq.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la  -ltrack -lapex -lmagcoord -lpomme -lcommon -lm -lmagfield -lmsynth -lbspline2 -lgsl -lnetcdf ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

print_fillion2_SOURCES = print_fillion2.c
print_fillion2_CFLAGS = -fopenmp
print_fillion2_LDFLAGS = -fopenmp
print_fillion2_LDADD = libpca3d.la $(top_builddir)/green/libgreen.la  $(top_builddir)/mageq/libmageq.la $(top_builddir)/lapack_wrapper/liblapack_wrapper.la  -ltrack -lapex -lmagcoord -lpomme -lcommon -lm -lmagfield -lmsynth -lbspline2 -lgsl -lnetcdf ${common_libs} -lshtns_omp -lfftw3 -lfftw3_omp

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

print_var3b_SOURCES = print_var3b.c
print_var3b_CFLAGS = -fopenmp
print_var3b_LDFLAGS = -fopenmp
print_var3b_LDADD = libpca3d.la -lcommon -lm -lmagfield -lgsl -lshtns_omp -lfftw3 -lfftw3_omp -lnetcdf ${common_libs}

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

#stage1c_SOURCES = stage1c.c
#stage1c_CFLAGS = -fopenmp
#stage1c_LDFLAGS = -fopenmp
#stage1c_LDADD = libpca3d.la -lcommon -lm -lmagfield ~/usr/lib/libgsl.a -lptcblas -latlas -lpthread -lshtns_omp -lfftw3 -lfftw3_omp -lnetcdf

print_SH_SOURCES = print_SH.c
print_SH_CFLAGS = -fopenmp
print_SH_LDFLAGS = -fopenmp
print_SH_LDADD = libpca3d.la -lcommon -lm -lmagfield ~/usr/lib/libgsl.a -lptcblas -latlas -lpthread -lshtns_omp -lfftw3 -lfftw3_omp -lnetcdf

