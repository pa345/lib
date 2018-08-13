#!/usr/bin/env gnuplot
#
# Plot time series of Stage2 residuals

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

# line types / colors
LT1 = 5
LT2 = 6

set out 'stage2.png'
set xrange ["1199145600":"1483228799"]

set multiplot layout nrow,ncol

set yrange [-200:200]
set key bottom right horizontal tc variable font "Helvetica Bold,18"

set title 'DMSP F-16 scalar residuals with CHAOS'
plot 'F16_stage1.txt' us 1:(abs($7) < 55 ? $12-$16 : 1/0) every nevery lt LT1 w dot ti "Stage1 F", \
     'F16_stage2.txt' us 1:(abs($7) < 55 ? $12-$16 : 1/0) every nevery lt LT2 w dot ti "Stage2 F", \

load 'incrow.cfg'

set title 'DMSP F-17 scalar residuals with CHAOS'
plot 'F17_stage1.txt' us 1:(abs($7) < 55 ? $12-$16 : 1/0) every nevery lt LT1 w dot ti "", \
     'F17_stage2.txt' us 1:(abs($7) < 55 ? $12-$16 : 1/0) every nevery lt LT2 w dot ti "", \

load 'incrow.cfg'

set title 'DMSP F-18 scalar residuals with CHAOS'
plot 'F18_stage1.txt' us 1:(abs($7) < 55 ? $12-$16 : 1/0) every nevery lt LT1 w dot ti "", \
     'F18_stage2.txt' us 1:(abs($7) < 55 ? $12-$16 : 1/0) every nevery lt LT2 w dot ti "", \

load 'inccolumn.cfg'

set title 'DMSP F-16 B_z residuals with CHAOS'
plot 'F16_stage1.txt' us 1:(abs($7) < 55 ? $11-$15 : 1/0) every nevery lt LT1 w dot ti "", \
     'F16_stage2.txt' us 1:(abs($7) < 55 ? $11-$15 : 1/0) every nevery lt LT2 w dot ti ""

load 'incrow.cfg'

set title 'DMSP F-17 B_z residuals with CHAOS'
plot 'F17_stage1.txt' us 1:(abs($7) < 55 ? $11-$15 : 1/0) every nevery lt LT1 w dot ti "", \
     'F17_stage2.txt' us 1:(abs($7) < 55 ? $11-$15 : 1/0) every nevery lt LT2 w dot ti ""

load 'incrow.cfg'

set title 'DMSP F-18 B_z residuals with CHAOS'
plot 'F18_stage1.txt' us 1:(abs($7) < 55 ? $11-$15 : 1/0) every nevery lt LT1 w dot ti "", \
     'F18_stage2.txt' us 1:(abs($7) < 55 ? $11-$15 : 1/0) every nevery lt LT2 w dot ti ""

unset multiplot
