bin_PROGRAMS = print print_lt print_plp calc_euler

commonlibs = -lmainlib -lbspline2 -lsofa_c -lcspice -lmwmclmcrrt -lmlchaos ~/usr/lib/libgsl.a -lgslcblas -lm -lcdf -lgfortran -lfftw3 -lz -lconfig ~/usr/lib/liblapacke.a ~/usr/lib/liblapack.a -lf77blas -lcblas -latlas

print_SOURCES = print.c
print_LDADD = ${commonlibs}

print_lt_SOURCES = print_lt.c
print_lt_LDADD = ${commonlibs}

print_plp_SOURCES = print_plp.c
print_plp_LDADD = ${commonlibs}

calc_euler_SOURCES = calc_euler.c
calc_euler_LDADD = ${commonlibs}
