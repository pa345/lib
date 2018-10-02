#!/usr/bin/env gnuplot
#
# Plot scalar rms values vs longitude

set term pngcairo enh col size 1600,1000 font ",18"
set out "rms_lon.png"

file = 'satrms.dat'
set yrange [0:100]
set point 1

set xlabel "longitude (degrees)"
set ylabel "scalar rms (nT)"
set title "F-17 scalar rms values for 2015"
plot file us 3:10 ti "Scalar RMS", 50 ti "Threshold"
