#!/usr/bin/env gnuplot
#
# Plot Euler angles

nrow = 3
ncol = 1

load 'multi_default.cfg'

plotheight = 2.0
plotwidth = 10.0
vbuffer = 0.8
l = 1.3

load 'multi_defs.cfg'
load 'multi_png.cfg'

satnum = '3'
iternum = '5'
outfile = 'euler.'.satnum.'_iter'.iternum.'.png'
set out outfile

str = sprintf('Generating %s', outfile)
print str

unset xlabel
unset key
load 'grid.cfg'
load 'lines2.cfg'
set xrange [2009:*]
mylw = 4

set multiplot layout nrow,ncol

file = '../output/euler.'.satnum.'.iter'.iternum.'.txt'

set ylabel "alpha (degrees)"
plot file us 2:3 w lp lw mylw lt 3

load 'incrow.cfg'

set ylabel "beta (degrees)"
plot file us 2:4 w lp lw mylw lt 3

load 'incrow.cfg'

set xlabel "time"
set ylabel "gamma (degrees)"
plot file us 2:5 w lp lw mylw lt 3

unset multiplot
