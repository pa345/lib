#!/usr/bin/env gnuplot
#
# Plot time series VFM residuals before/after spike correction

nrow=3
ncol=1

load 'multi_default.cfg'

plotwidth = 10.0
fontsize = ",12"
vbuffer = 0.8

load 'multi_defs.cfg'
load 'multi_png.cfg'

set out "spikes.png"

set ylabel "nT"
load 'xtimeon.cfg'

file = 'stage2_spikes.dat'

set key tc variable
idx1 = 8
idx2 = 11

set multiplot layout nrow,ncol

set yrange [-150:150]
set title "VFM X spike correction"
plot file us 1:4 index idx1:idx2 ti "Original VFM X", \
     file us 1:14 index idx1:idx2 ti "Upper limit", \
     file us 1:15 index idx1:idx2 ti "Lower limit", \
     file us 1:23 index idx1:idx2 ti "Corrected VFM X", \
     file us 1:($20 == 1 ? $4 : 1/0) index idx1:idx2 pt 4 ps 1 ti "Detected outliers"

load 'incrow.cfg'

set title "VFM Y spike correction"
plot file us 1:5 index idx1:idx2 ti "Original VFM Y", \
     file us 1:16 index idx1:idx2 ti "Upper limit", \
     file us 1:17 index idx1:idx2 ti "Lower limit", \
     file us 1:24 index idx1:idx2 ti "Corrected VFM Y", \
     file us 1:($21 == 1 ? $5 : 1/0) index idx1:idx2 pt 4 ps 1 ti "Detected outliers"

load 'incrow.cfg'

set title "VFM Z spike correction"
plot file us 1:6 index idx1:idx2 ti "Original VFM Z", \
     file us 1:18 index idx1:idx2 ti "Upper limit", \
     file us 1:19 index idx1:idx2 ti "Lower limit", \
     file us 1:25 index idx1:idx2 ti "Corrected VFM Z", \
     file us 1:($22 == 1 ? $6 : 1/0) index idx1:idx2 pt 4 ps 1 ti "Detected outliers"

unset multiplot
