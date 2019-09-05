check_PROGRAMS = main

common_libs = -lmainlib -lbspline2 -lsofa_c -lcspice -lmwmclmcrrt -lmlchaos ~/usr/lib/libgsl.a -lm -lcdf -lgfortran -lfftw3 -lz -lconfig ~/usr/lib/liblapacke.a ~/usr/lib/liblapack.a -lptf77blas -lptcblas -latlas -lgfortran -lpthread

main_SOURCES = main.c
main_LDADD = ${common_libs}
