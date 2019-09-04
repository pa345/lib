bin_PROGRAMS = print print_all print_efi print_lp print_lt print_crust print_alt print_tracks residuals plotsep swarmcdf2fday

commonlibs = -lmainlib -lbspline2 -lsofa_c -lcspice -lmwmclmcrrt -lmlchaos ~/usr/lib/libgsl.a -lm -lcdf -lgfortran -lfftw3 -lz -lconfig ~/usr/lib/liblapack.a ~/usr/lib/libf77blas.a ~/usr/lib/libcblas.a ~/usr/lib/libatlas.a

plotsep_SOURCES = plotsep.c
plotsep_LDADD = ${commonlibs}

print_SOURCES = print.c
print_LDADD = ${commonlibs}

print_tracks_SOURCES = print_tracks.c
print_tracks_LDADD = ${commonlibs}

print_lt_SOURCES = print_lt.c
print_lt_LDADD = ${commonlibs}

print_all_SOURCES = print_all.c
print_all_LDADD = ${commonlibs}

print_alt_SOURCES = print_alt.c
print_alt_LDADD = ${commonlibs}

print_efi_SOURCES = print_efi.c
print_efi_LDADD = ${commonlibs}

print_lp_SOURCES = print_lp.c
print_lp_LDADD = ${commonlibs}

print_crust_SOURCES = print_crust.c
print_crust_LDADD = ${commonlibs}

residuals_SOURCES = residuals.c
residuals_LDFLAGS = -fopenmp
residuals_LDADD = ${commonlibs}

swarmcdf2fday_SOURCES = swarmcdf2fday.c
swarmcdf2fday_LDADD = ${commonlibs}

AM_CPPFLAGS = -I$(top_builddir)/pca -I$(top_builddir)/green
