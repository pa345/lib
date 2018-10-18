bin_PROGRAMS = print print_all print_efi print_lp print_lt print_crust print_alt print_tracks residuals plotsep swarmcdf2fday

commonlibs = $(top_builddir)/estist/libestist_calc.la -lm ~/usr/lib/liblapack.a ~/usr/lib/libf77blas.a ~/usr/lib/libatlas.a -ltrack -lpomme -lapex -lcommon -lmsynth -lbspline2 -lcdf -lsatdata -lindices ~/usr/lib/libgsl.a -lgslcblas -lgfortran

plotsep_SOURCES = plotsep.c
plotsep_LDADD = -lapex -lcommon -lmsynth -lbspline2 -lm -lgfortran ${commonlibs}

print_SOURCES = print.c
print_LDADD = -lm ${commonlibs}

print_tracks_SOURCES = print_tracks.c
print_tracks_LDADD = -lm ${commonlibs}

print_lt_SOURCES = print_lt.c
print_lt_LDADD = -lm ${commonlibs}

print_all_SOURCES = print_all.c
print_all_LDADD = -lapex $(top_builddir)/euler/libeuler.la -lm ${commonlibs}

print_alt_SOURCES = print_alt.c
print_alt_LDADD = -lapex -lm ${commonlibs}

print_efi_SOURCES = print_efi.c
print_efi_LDADD = -lapex -lm ${commonlibs} -lgfortran

print_lp_SOURCES = print_lp.c
print_lp_LDADD = -lapex -lm ${commonlibs} -lgfortran

print_crust_SOURCES = print_crust.c
print_crust_LDADD = -lapex -lm ${commonlibs}

residuals_SOURCES = residuals.c
residuals_LDFLAGS = -fopenmp
residuals_LDADD = ${commonlibs}

swarmcdf2fday_SOURCES = swarmcdf2fday.c
swarmcdf2fday_LDADD = -lm ${commonlibs}

AM_CPPFLAGS = -I$(top_builddir)/euler -I$(top_builddir)/chaos -I$(top_builddir)/magdata -I$(top_builddir)/pca -I$(top_builddir)/green -I$(top_builddir)/estist
