#!/usr/bin/env gnuplot

nrow = 2
ncol = 2

load 'multi_default.cfg'

plotwidth = "5.0"
b = 0.3
hbuffer = 1.8
r = -0.7

load 'multi_defs.cfg'
load 'multi_png.cfg'

set out "B_alt.png"

file = 'data_alt.txt'

set pm3d map interp 0,0
set xrange [-25:25]
set yrange [*:300]
set xtics -20,5,30
set ylabel "altitude (km)"
unset key

set multiplot layout nrow,ncol

set cblabel "nT"

set format x ""
set title "B_r (longitude = 150 degrees)"
splot file us 1:2:6

load 'incrow.cfg'

set format x "%g"
set xlabel "geographic latitude (degrees)"

set title "B_{/Symbol \152}"
splot file us 1:2:8

unset xlabel
set format x ""

load 'inccolumn.cfg'

set title "B_{/Symbol \161}"
splot file us 1:2:7

load 'incrow.cfg'

set format x "%g"
set xlabel "geographic latitude (degrees)"

set title "|B|"
splot file us 1:2:9


unset multiplot
