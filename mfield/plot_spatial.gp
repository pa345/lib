#!/home/palken/usr/bin/gnuplot
#
# Plot time/latitude scatter plot of available data

nrow = 4
ncol = 1

load 'multi_default.cfg'

plotwidth = 16.0
b = 0.2
hbuffer = 0.0

load 'multi_defs.cfg'
load 'multi_png.cfg'

nskip = 3

set style fill transparent solid 1.0 noborder
set style circle radius 0.004

load 'lines2.cfg'
load 'ylaton.cfg'
load 'grid.cfg'
#set xrange [2013.5:*]

outfile = 'output/spatial.png'
set out outfile

set key top left inside tc variable

set multiplot layout nrow,ncol

set format x ""

set title "X component"
plot 'output/datamap0_X.dat' us 1:3 every nskip w p ti "Swarm A", \
     'output/datamap1_X.dat' us 1:3 every nskip w p lt 3 ti "Swarm B"

load 'incrow.cfg'

set title "Y component"
plot 'output/datamap0_Y.dat' us 1:3 every nskip w p ti "", \
     'output/datamap1_Y.dat' us 1:3 every nskip w p lt 3 ti ""

load 'incrow.cfg'

set title "Z component"
plot 'output/datamap0_Z.dat' us 1:3 every nskip w p ti "", \
     'output/datamap1_Z.dat' us 1:3 every nskip w p lt 3 ti ""

load 'incrow.cfg'

set format x "%g"
set xlabel "time"

set title "F component"
plot 'output/datamap0_F.dat' us 1:3 every nskip w p ti "", \
     'output/datamap1_F.dat' us 1:3 every nskip w p lt 3 ti ""

set format x ""
unset xlabel

unset multiplot

outstr = sprintf('output is %s...', outfile)
print outstr

