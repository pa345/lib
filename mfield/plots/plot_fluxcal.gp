#!/usr/bin/env gnuplot
#
# Plot fluxgate calibration parameters

nrow = 3
ncol = 1

load 'multi_default.cfg'

plotheight = 2.0
plotwidth = 10.0
vbuffer = 0.8

load 'multi_defs.cfg'
load 'multi_png.cfg'

satnum = '3'
iternum = '5'
outfile = 'fluxcal.'.satnum.'_iter'.iternum.'.png'
set out outfile

str = sprintf('Generating %s', outfile)
print str

unset xlabel
load 'grid.cfg'
load 'lines2.cfg'

set multiplot layout nrow,ncol

file = '../output/fluxcal.'.satnum.'.iter'.iternum.'.txt'
set title "Scale factors"
plot file us 2:3 w lp ti "S_x", \
     file us 2:4 w lp ti "S_y", \
     file us 2:5 w lp ti "S_z"

load 'incrow.cfg'

set ylabel "nT"
set title "Offsets"
plot file us 2:6 w lp ti "O_x", \
     file us 2:7 w lp ti "O_y", \
     file us 2:8 w lp ti "O_z"

load 'incrow.cfg'

set xlabel "time"
set ylabel "degrees"
set title "Non-orthogonality angles"
plot file us 2:9 w lp ti "U_1", \
     file us 2:10 w lp ti "U_2", \
     file us 2:11 w lp ti "U_3"

unset multiplot