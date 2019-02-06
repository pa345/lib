#!/usr/bin/env gnuplot

nrow = 4
ncol = 3

load 'multi_default.cfg'

plotwidth = 8.0
hbuffer = 1.0
b = 0.5

load 'multi_defs.cfg'

#load 'multi_png.cfg'
#set out "EZIE2.png"
load 'multi_pdf.cfg'
set out "EZIE2.pdf"

file_signal = 'signal.dat'
file_fitted = 'fitted.dat'
file_fitted_exact = 'fitted.dat.0'

load 'xyborder.cfg'
load 'grid.cfg'
set ylabel "nT"
set key top left

set multiplot layout nrow,ncol

set format x ""
set yrange [-400:600]

plot file_signal index 0 us 1:2 w li lw 4 lc rgb "blue" ti "Exact B_r", \
     file_signal index 0 us 1:3 w li lw 2 lc rgb "blue" ti "Noisy B_r", \
     file_fitted index 0 us 1:2 w li lw 4 dt 2 lc rgb "blue" ti "Fitted B_r"

load 'incrow.cfg'

plot file_signal index 1 us 1:2 w li lw 4 lc rgb "orange" ti "", \
     file_signal index 1 us 1:3 w li lw 2 lc rgb "orange" ti "", \
     file_fitted index 1 us 1:2 w li lw 4 dt 2 lc rgb "orange" ti ""

load 'incrow.cfg'

plot file_signal index 2 us 1:2 w li lw 4 lc rgb "dark-green" ti "", \
     file_signal index 2 us 1:3 w li lw 2 lc rgb "dark-green" ti "", \
     file_fitted index 2 us 1:2 w li lw 4 dt 2 lc rgb "dark-green" ti ""

load 'incrow.cfg'
set format x "%g"
set xlabel "latitude (degrees)"

plot file_signal index 3 us 1:2 w li lw 4 lc rgb "red" ti "", \
     file_signal index 3 us 1:3 w li lw 2 lc rgb "red" ti "", \
     file_fitted index 3 us 1:2 w li lw 4 dt 2 lc rgb "red" ti ""

load 'inccolumn.cfg'

set format x ""
unset xlabel
set yrange [-400:100]

plot file_signal index 0 us 1:4 w li lw 4 lc rgb "blue" ti "Exact B_t", \
     file_signal index 0 us 1:5 w li lw 2 lc rgb "blue" ti "Noisy B_t", \
     file_fitted index 0 us 1:3 w li lw 4 dt 2 lc rgb "blue" ti "Fitted B_t"

load 'incrow.cfg'

plot file_signal index 1 us 1:4 w li lw 4 lc rgb "orange" ti "", \
     file_signal index 1 us 1:5 w li lw 2 lc rgb "orange" ti "", \
     file_fitted index 1 us 1:3 w li lw 4 dt 2 lc rgb "orange" ti ""

load 'incrow.cfg'

plot file_signal index 2 us 1:4 w li lw 4 lc rgb "dark-green" ti "", \
     file_signal index 2 us 1:5 w li lw 2 lc rgb "dark-green" ti "", \
     file_fitted index 2 us 1:3 w li lw 4 dt 2 lc rgb "dark-green" ti ""

load 'incrow.cfg'
set format x "%g"
set xlabel "latitude (degrees)"

plot file_signal index 3 us 1:4 w li lw 4 lc rgb "red" ti "", \
     file_signal index 3 us 1:5 w li lw 2 lc rgb "red" ti "", \
     file_fitted index 3 us 1:3 w li lw 4 dt 2 lc rgb "red" ti ""

load 'inccolumn.cfg'

set format x ""
unset xlabel
set ylabel "mA/m"
set yrange [-50:500]

plot file_fitted_exact index 0 us 1:4 w li lw 4 lc rgb "blue" ti "Exact J_p", \
     file_fitted index 0 us 1:($4-0) w li lw 4 dt 2 lc rgb "blue" ti "Fitted J_p"

load 'incrow.cfg'

plot file_fitted_exact index 1 us 1:4 w li lw 4 lc rgb "orange" ti "", \
     file_fitted index 1 us 1:($4-0) w li lw 4 dt 2 lc rgb "orange" ti ""

load 'incrow.cfg'

plot file_fitted_exact index 2 us 1:4 w li lw 4 lc rgb "dark-green" ti "", \
     file_fitted index 2 us 1:($4-0) w li lw 4 dt 2 lc rgb "dark-green" ti ""

load 'incrow.cfg'
set format x "%g"
set xlabel "latitude (degrees)"

plot file_fitted_exact index 3 us 1:4 w li lw 4 lc rgb "red" ti "", \
     file_fitted index 3 us 1:($4-0) w li lw 4 dt 2 lc rgb "red" ti ""

unset multiplot
