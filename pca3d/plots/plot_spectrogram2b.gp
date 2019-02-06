#!/usr/bin/env gnuplot
#
# Plot Jr/Jt/Jp time series along with spectrograms

nrow = 3
ncol = 1

load 'multi_default.cfg'

plotheight = 1.5
plotwidth = 5.5
fontsize = ",10"
r = -0.9
l = 0.5
hbuffer = 1.3

load 'multi_defs.cfg'
load 'multi_png.cfg'

set output "spectrogram2b.png"

file_freq = 'spectrogram2b.txt'

load 'lines2.cfg'
load 'xtimeon.cfg'
load 'jet.pal'

# xtics every ndays days
ndays = 5
set xtics 3600*24*ndays

set palette maxcol 0

set multiplot layout nrow,ncol

unset key
set cblabel "Power ({/Symbol \155}A/m^2)"
set ylabel "Period (days)"

set pm3d map interp 0,0

set title "Spectrogram of qtlm"
splot file_freq us 1:3:4

load 'incrow.cfg'

set title "Spectrogram of plm"
splot file_freq us 1:3:5

load 'incrow.cfg'

set title "Spectrogram of qlm"
splot file_freq us 1:3:6

unset multiplot
