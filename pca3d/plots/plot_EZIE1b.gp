#!/usr/bin/env gnuplot

#set term pngcairo enh col
#set out "EZIE1b.png"
set term pdfcairo enh col
set out "EZIE1b.pdf"

file = 'data_map.txt'

set pm3d map interp 0,0

set xrange [138:156]
set yrange [-10:20]
unset key

set xlabel "longitude (degrees)"
set ylabel "latitude (degrees)"
set cblabel "nT"

load 'jet.pal'
set palette maxcol 0

set size 0.96

# factor to account for TIEGCM not provided below 109km
fac = 2.2

mylw = 4
#lon1 = 136.5
#lon2 = 141.0
#lon3 = 142.5
#lon4 = 145.5

lon1 = 143.5
lon2 = 148.5
lon3 = 150.0
lon4 = 152.5

set arrow 1 from lon1,-10 to lon1,20 nohead front lw mylw dt "-" lc rgb "blue"
set arrow 2 from lon2,-10 to lon2,20 nohead front lw mylw dt "-" lc rgb "orange"
set arrow 3 from lon3,-10 to lon3,20 nohead front lw mylw dt "-" lc rgb "dark-green"
set arrow 4 from lon4,-10 to lon4,20 nohead front lw mylw dt "-" lc rgb "red"

set title "B_{/Symbol \161}"
#set title "B_r"
splot file us 1:2:($7*fac)