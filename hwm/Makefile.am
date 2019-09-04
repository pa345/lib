lib_LTLIBRARIES = libhwm.la

libhwm_la_SOURCES = hwm.c hwm14.f90

AM_CPPFLAGS =

FC = gfortran

AM_FCFLAGS = -g -O2
AM_CFLAGS = -std=c99 -g -O2 -Wall -W -D_XOPEN_SOURCE

check_PROGRAMS = test checkhwm14

test_SOURCES = test.c
test_LDADD = libhwm.la -lm -lgfortran -lindices

checkhwm14_SOURCES = checkhwm14.f90
checkhwm14_LDADD = libhwm.la -lm -lgfortran -lindices
