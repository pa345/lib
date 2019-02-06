#!/usr/bin/env gnuplot
#
# Plot PCA modes from SH analysis

nrow = 3
ncol = 2

load 'multi_default.cfg'

plotheight = 2.0
plotwidth = 3.0
fontsize = ",10"
r = -0.3
l = 0.6
hbuffer = 0.8

load 'multi_defs.cfg'
load 'multi_png.cfg'

set out "modes.png"

unset key
set pm3d map interp 0,0
set palette maxcol 0
load 'jet.pal'

set xrange [-180:180]
set xtics -180,45,180
load 'ylaton.cfg'

set multiplot layout nrow,ncol

set title "J_{/Symbol \161}, Mode 1"
splot 'J_01.txt' us 1:2:4

load 'incrow.cfg'

set title "J_{/Symbol \161}, Mode 2"
splot 'J_02.txt' us 1:2:4

load 'incrow.cfg'

set xlabel "longitude (degrees)"

set title "J_{/Symbol \161}, Mode 3"
splot 'J_03.txt' us 1:2:4

unset xlabel

load 'inccolumn.cfg'

set title "J_{/Symbol \152}, Mode 1"
splot 'J_01.txt' us 1:2:5

load 'incrow.cfg'

set title "J_{/Symbol \152}, Mode 2"
splot 'J_02.txt' us 1:2:5

load 'incrow.cfg'

set xlabel "longitude (degrees)"

set title "J_{/Symbol \152}, Mode 3"
splot 'J_03.txt' us 1:2:5

unset multiplot
