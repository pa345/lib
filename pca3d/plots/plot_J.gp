#!/usr/bin/env gnuplot

nrow = 2
ncol = 2

load 'multi_default.cfg'

plotwidth = "5.0"
b = 0.3
hbuffer = 2.0
r = -0.7

load 'multi_defs.cfg'
load 'multi_png.cfg'

set out "J_alt.png"

file = 'data_alt.txt'

set pm3d map interp 0,0
set xrange [-25:25]
set yrange [*:300]
set xtics -20,5,30
set ylabel "altitude (km)"
unset key

set multiplot layout nrow,ncol

set cblabel "uA/m^2"

set format x ""
set title "J_r (longitude = 150 degrees)"
splot file us 1:2:($3*1e6)

load 'incrow.cfg'

set format x "%g"
set xlabel "geographic latitude (degrees)"

set title "J_{/Symbol \152}"
splot file us 1:2:($5*1e6)

unset xlabel
set format x ""

load 'inccolumn.cfg'

set format x "%g"
set xlabel "geographic latitude (degrees)"

set title "J_{/Symbol \161}"
splot file us 1:2:($4*1e6)

unset multiplot
