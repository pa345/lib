#!/usr/bin/env gnuplot
#
# Plot PCA modes from SH analysis

nrow = 3
ncol = 1

load 'multi_default.cfg'

plotheight = 2.0
plotwidth = 3.0
fontsize = ",10"
r = -0.3
l = 0.6
hbuffer = 0.8

load 'multi_defs.cfg'
load 'multi_png.cfg'

mode = ARG1

outfile = 'mode_'.mode.'.png'
set out outfile
file = "J_".mode.".txt"

unset key
set pm3d map interp 0,0
set palette maxcol 0
load 'jet.pal'

set xrange [-180:180]
set xtics -180,45,180
load 'ylaton.cfg'

set multiplot layout nrow,ncol

set title "J_r, Mode ".mode
splot file us 1:2:3

load 'incrow.cfg'

set title "J_{/Symbol \161}, Mode ".mode
splot file us 1:2:4

load 'incrow.cfg'

set xlabel "longitude (degrees)"

set title "J_{/Symbol \152}, Mode ".mode
splot file us 1:2:5

unset multiplot

str = sprintf('output is %s', outfile)
print str
