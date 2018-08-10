#!/usr/bin/env gnuplot
#
# Plot time series of Stage1 residuals

nrow=4
ncol=2

load 'multi_default.cfg'

plotwidth = 10.0
fontsize = ",12"
vbuffer = 0.7

load 'multi_defs.cfg'
load 'multi_png.cfg'

load 'xtimeon.cfg'
set ylabel "nT"

nevery = 3

file_f15 = 'F15_stage1.txt'
file_f16 = 'F16_stage1.txt'
file_f17 = 'F17_stage1.txt'
file_f18 = 'F18_stage1.txt'

# line types / colors
LT1 = 5
LT2 = 6

set out 'stage1a.png'
set xrange ["1199145600":"1483228799"]

set multiplot layout nrow,ncol

set yrange [-200:200]
set key bottom right horizontal tc variable font "Helvetica Bold,18"

set title 'DMSP F-15 scalar residuals with CHAOS'
plot file_f15 us 1:($12-$16) every nevery lt LT1 w dot ti "F"

load 'incrow.cfg'

set title 'DMSP F-16 scalar residuals with CHAOS'
plot file_f16 us 1:($12-$16) every nevery lt LT1 w dot ti ""

load 'incrow.cfg'

set title 'DMSP F-17 scalar residuals with CHAOS'
plot file_f17 us 1:($12-$16) every nevery lt LT1 w dot ti ""

load 'incrow.cfg'

set title 'DMSP F-18 scalar residuals with CHAOS'
plot file_f18 us 1:($12-$16) every nevery lt LT1 w dot ti ""

load 'inccolumn.cfg'

set title 'DMSP F-15 B_z residuals with CHAOS'
plot file_f15 us 1:($11-$15) every nevery lt LT1 w dot ti "All B_z"

load 'incrow.cfg'

set title 'DMSP F-16 B_z residuals with CHAOS'
plot file_f16 us 1:($11-$15) every nevery lt LT1 w dot ti ""

load 'incrow.cfg'

set title 'DMSP F-17 B_z residuals with CHAOS'
plot file_f17 us 1:($11-$15) every nevery lt LT1 w dot ti ""

load 'incrow.cfg'

set title 'DMSP F-18 B_z residuals with CHAOS'
plot file_f18 us 1:($11-$15) every nevery lt LT1 w dot ti ""

unset multiplot
