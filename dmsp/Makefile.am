bin_PROGRAMS = print print_LTAN stage0 stage2a stage2b atest

commonlibs = -lmainlib -lbspline2 -lsofa_c -lcspice -lmwmclmcrrt -lmlchaos ~/usr/lib/libgsl.a -lgslcblas -lm -lcdf -lgfortran -lfftw3 -lz -lconfig ~/usr/lib/liblapacke.a ~/usr/lib/liblapack.a -lf77blas -latlas

AM_CPPFLAGS =

stage0_SOURCES = stage0.c
stage0_LDADD = ${commonlibs}

atest_SOURCES = atest.c
atest_LDADD = ${commonlibs}

stage2a_SOURCES = jump.c fluxcal.c magcal.c stage2a.c
stage2a_LDFLAGS = -fopenmp
stage2a_LDADD = ${commonlibs}

stage2b_SOURCES = stage2b.c magcal.c
stage2b_LDADD = ${commonlibs}

print_SOURCES = print.c
print_LDADD = ${commonlibs}

print_LTAN_SOURCES = print_LTAN.c
print_LTAN_LDADD = ${commonlibs}
