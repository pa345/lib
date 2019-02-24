#!/usr/bin/env gnuplot
#
# Plot time series (histogram) of model data distribution
# Use plothist.sh to generate data files

### n: change this parameter to equal the number of data sets to be plotted
nrow = 3
ncol = 2

load 'multi_default.cfg'

hbuffer = 1.2
r = 0.0

load 'multi_defs.cfg'
load 'multi_png.cfg'
set out "hist.png"

load 'lines2.cfg'

set ylabel "track scalar rms (nT)"
#set xrange [2013.5:*]
set yrange [0:*]
set ytics nomirror

set multiplot layout nrow,ncol

set xtics 0.5 format ""

set ylabel "vector data per 10 days"

set title "Temporal distribution of X vector data for main field modeling"
plot 'hist_vec_X.sat0' us 1:3 w boxes fill so ti "Swarm A", \
     'hist_vec_X.sat1' us 1:3 w boxes fill so ti "Swarm B"

unset key

load 'incrow.cfg'

set title "Temporal distribution of Y vector data for main field modeling"
plot 'hist_vec_Y.sat0' us 1:3 w boxes fill so, \
     'hist_vec_Y.sat1' us 1:3 w boxes fill so

load 'incrow.cfg'

set xtics format "%g"

set title "Temporal distribution of Z vector data for main field modeling"
plot 'hist_vec_Z.sat0' us 1:3 w boxes fill so, \
     'hist_vec_Z.sat1' us 1:3 w boxes fill so

set xtics 0.5 format ""

load 'inccolumn.cfg'

set title "Temporal distribution of scalar data for main field modeling"
set ylabel "scalar data per 10 days"
plot 'hist_scal.sat0' us 1:3 w boxes fill so ti "Swarm A", \
     'hist_scal.sat1' us 1:3 w boxes fill so ti "Swarm B"

load 'incrow.cfg'

set xtics format "%g"

set title "Temporal distribution of vector data for Euler angles"
set ylabel "Euler angle data per 10 days"
plot 'hist_euler.sat0' us 1:3 w boxes fill so, \
     'hist_euler.sat1' us 1:3 w boxes fill so

unset multiplot
