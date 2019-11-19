#!/usr/bin/env gnuplot
#

nrow = 4
ncol = 1

load 'multi_default.cfg'

plotheight = 2.0
plotwidth = 10.0
vbuffer = 0.8

load 'multi_defs.cfg'
load 'multi_png.cfg'

filedir = '../output'
satnum = '1'
iternum = '15'
set out 'res'.satnum.'_iter'.iternum.'.png'

unset key
nskip = 1

load 'xlaton.cfg'

unset xlabel
set ylabel "nT"
set yrange [-50:50]
set ytics 25
#set yrange [-400:400]

set multiplot layout nrow,ncol

file = filedir.'/res'.satnum.'_X_iter'.iternum.'.dat'
set title "X residuals"
plot file us 5:13 every nskip w p

load 'incrow.cfg'

file = filedir.'/res'.satnum.'_Y_iter'.iternum.'.dat'
set title "Y residuals"
plot file us 5:13 every nskip w p

load 'incrow.cfg'

file = filedir.'/res'.satnum.'_Z_iter'.iternum.'.dat'
set title "Z residuals"
plot file us 5:13 every nskip w p

load 'incrow.cfg'

set xlabel "QD latitude (degrees)"

file = filedir.'/res'.satnum.'_F_iter'.iternum.'.dat'
set title "F residuals"
plot file us 5:12 every nskip w p

unset multiplot
