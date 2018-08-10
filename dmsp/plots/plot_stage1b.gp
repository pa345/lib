#!/usr/bin/env gnuplot
#
# Plot Stage1 residuals vs QD latitude

nrow=6
ncol=2

load 'multi_default.cfg'

plotwidth = 10.0
fontsize = ",12"
vbuffer = 0.7

load 'multi_defs.cfg'
load 'multi_png.cfg'

load 'xlaton.cfg'
unset xlabel

set ylabel "nT"

nevery = 1

# line types / colors
LT1 = 5
LT2 = 6

set key bottom right horizontal tc variable font "Helvetica Bold,18"

do for [satidx=15:15] {

  file = 'F'.satidx.'_stage1.txt'
  outfile = 'stage1_F'.satidx.'.png'
  set out outfile

  str = sprintf('Generating %s...', outfile)
  print str

  set yrange [-150:150]

  set multiplot layout nrow,ncol

  do for [year=2008:2013] {

    tmin = year
    tmax = year + 1

    if (year > 2008) { unset key }
    if (year == 2013) { set xlabel "QD latitude (degrees)" }

    set title 'DMSP F-15 scalar residuals with CHAOS between ['.tmin.','.tmax.']'
    plot file us 7:(($2 > tmin & $2 < tmax & $8 == 1) ? $12-$16 : 1/0) every nevery lt LT1 w dot ti "Ascending F", \
         file us 7:(($2 > tmin & $2 < tmax & $8 == -1) ? $12-$16 : 1/0) every nevery lt LT2 w dot ti "Descending F"

    load 'incrow.cfg'
  }

  load 'inccolumn.cfg'
  set key
  unset xlabel

  do for [year=2008:2013] {

    tmin = year
    tmax = year + 1

    if (year > 2008) { unset key }
    if (year == 2013) { set xlabel "QD latitude (degrees)" }

    set title 'DMSP F-15 B_z residuals with CHAOS between ['.tmin.','.tmax.']'
    plot file us 7:(($2 > tmin & $2 < tmax & $8 == 1) ? $11-$15 : 1/0) every nevery lt LT1 w dot ti "Ascending B_z", \
         file us 7:(($2 > tmin & $2 < tmax & $8 == -1) ? $11-$15 : 1/0) every nevery lt LT2 w dot ti "Descending B_z"

    load 'incrow.cfg'
  }

  unset multiplot

  load 'multi_reset.cfg'
}

do for [satidx=16:18] {

  file = 'F'.satidx.'_stage1.txt'
  outfile = 'stage1_F'.satidx.'.png'
  set out outfile

  str = sprintf('Generating %s...', outfile)
  print str

  set yrange [-150:150]

  set multiplot layout nrow,ncol

  do for [year=2009:2016] {

    tmin = year
    tmax = year + 1

    if (year > 2009) { unset key }
    if (year == 2016) { set xlabel "QD latitude (degrees)" }

    set title 'DMSP F-'.satidx.' scalar residuals with CHAOS between ['.tmin.','.tmax.']'
    plot file us 7:(($2 > tmin & $2 < tmax & $8 == 1) ? $12-$16 : 1/0) every nevery lt LT1 w dot ti "Ascending F", \
         file us 7:(($2 > tmin & $2 < tmax & $8 == -1) ? $12-$16 : 1/0) every nevery lt LT2 w dot ti "Descending F"

    load 'incrow.cfg'
  }

  load 'inccolumn.cfg'
  set key
  unset xlabel

  do for [year=2009:2016] {

    tmin = year
    tmax = year + 1

    if (year > 2009) { unset key }
    if (year == 2016) { set xlabel "QD latitude (degrees)" }

    set title 'DMSP F-'.satidx.' B_z residuals with CHAOS between ['.tmin.','.tmax.']'
    plot file us 7:(($2 > tmin & $2 < tmax & $8 == 1) ? $11-$15 : 1/0) every nevery lt LT1 w dot ti "Ascending B_z", \
         file us 7:(($2 > tmin & $2 < tmax & $8 == -1) ? $11-$15 : 1/0) every nevery lt LT2 w dot ti "Descending B_z"

    load 'incrow.cfg'
  }

  unset multiplot

  load 'multi_reset.cfg'
}
