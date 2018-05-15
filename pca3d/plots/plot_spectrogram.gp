#!/usr/bin/env gnuplot
#
# Plot Jr/Jt/Jp time series along with spectrograms

nrow = 3
ncol = 2

load 'multi_default.cfg'

plotheight = 1.5
plotwidth = 5.5
fontsize = ",10"
r = -0.9
l = 0.5
hbuffer = 1.3

load 'multi_defs.cfg'
load 'multi_png.cfg'

set output "spectrogram.png"

file_time = 'spectrogram_time.txt'
file_freq = 'spectrogram.txt'

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

set title "Spectrogram of J_r"
splot file_freq us 1:2:4

load 'incrow.cfg'

set title "Spectrogram of J_{/Symbol \161}"
splot file_freq us 1:2:5

load 'incrow.cfg'

set title "Spectrogram of J_{/Symbol \152}"
splot file_freq us 1:2:6

load 'inccolumn.cfg'

load 'grid.cfg'
unset pm3d
set ylabel "{/Symbol \155}A/m^2"

set title "J_r"
plot file_time us 1:2 w li

load 'incrow.cfg'

set title "J_{/Symbol \161}"
plot file_time us 1:3 w li

load 'incrow.cfg'

set title "J_{/Symbol \152}"
plot file_time us 1:4 w li

unset multiplot
