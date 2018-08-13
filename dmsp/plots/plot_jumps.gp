#!/usr/bin/env gnuplot
#
# Plot time series VFM residuals before/after jump correction

nrow=3
ncol=1

load 'multi_default.cfg'

plotwidth = 10.0
fontsize = ",12"
vbuffer = 0.8

load 'multi_defs.cfg'
load 'multi_png.cfg'

set out "jumps.png"

set ylabel "nT"
load 'xtimeon.cfg'

file = 'stage2_jumps.dat'

mylw = 4
idx1 = 8
idx2 = 11

set multiplot layout nrow,ncol

set yrange [-150:150]
set title "VFM X residuals"
plot file us 1:4 index idx1:idx2 w li lw mylw ti "Original VFM X", \
     file us 1:7 index idx1:idx2 w li lw mylw axes x1y2 ti "Test statistic", \
     file us 1:($10*10) index idx1:idx2 w li lw mylw axes x1y2 ti "Jump detected", \
     file us 1:13 index idx1:idx2 w li lw mylw lt 7 ti "Corrected VFM X"

load 'incrow.cfg'

set yrange [-100:100]
set title "VFM Y residuals"
plot file us 1:5 index idx1:idx2 w li lw mylw ti "Original VFM Y", \
     file us 1:8 index idx1:idx2 w li lw mylw axes x1y2 ti "Test statistic", \
     file us 1:($11*10) index idx1:idx2 w li lw mylw axes x1y2 ti "Jump detected", \
     file us 1:14 index idx1:idx2 w li lw mylw lt 7 ti "Corrected VFM Y"

load 'incrow.cfg'

set yrange [-100:100]
set title "VFM Z residuals"
plot file us 1:6 index idx1:idx2 w li lw mylw ti "Original VFM Z", \
     file us 1:9 index idx1:idx2 w li lw mylw axes x1y2 ti "Test statistic", \
     file us 1:($12*10) index idx1:idx2 w li lw mylw axes x1y2 ti "Jump detected", \
     file us 1:15 index idx1:idx2 w li lw mylw lt 7 ti "Corrected VFM Z"

unset multiplot
