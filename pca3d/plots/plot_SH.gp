#!/usr/bin/env gnuplot
#
# Plot SH time series

nrow = 2
ncol = 1

load 'multi_default.cfg'

plotheight = 1.5
plotwidth = 5.5
fontsize = ",10"
r = -0.1
l = 0.6

load 'multi_defs.cfg'
load 'multi_png.cfg'

set out "plot.png"

file = 'dat'

load 'grid.cfg'
load 'lines2.cfg'
load 'xtimeon.cfg'

set ylabel "nT"

set multiplot layout nrow,ncol

plot file us 1:2 w li ti "qt10"

load 'incrow.cfg'

plot file us 1:4 w li lt 2 ti "p10", \
     file us 1:6 w li lt 3 ti "q10"

unset multiplot
