#!/home/palken/usr/bin/gnuplot
#
# Plot F2 curves and corresponding L-curves for Sq fit
#
# Adjust 'logdir' and 'outdir' below

logdir = 'log_A'
outdir = 'log_A/plots'

cmd = sprintf('mkdir -p %s', outdir)
tstr = system(cmd)

fileF2 = logdir.'/F2.dat'
fileEEJ = logdir.'/EEJ.dat'
fileSqLcurve = logdir.'/Sq_lcurve.dat'
fileSqLcorner = logdir.'/Sq_lcorner.dat'
fileEEJLcurve = logdir.'/EEJ_lcurve.dat'
fileEEJLcorner = logdir.'/EEJ_lcorner.dat'
fileprof = logdir.'/profile.dat'

# Read timestamps into array
set term dumb
tarr=""
phiarr=""
ltarr=""
kparr=""
dirarr=""
SqR2arr=""
EEJR2arr=""
plot fileprof us ( tarr=tarr.stringcolumn(2).' ', \
                   phiarr=phiarr.stringcolumn(3).' ', \
                   ltarr=ltarr.stringcolumn(4).' ', \
                   kparr=kparr.stringcolumn(8).' ', \
                   dirarr=dirarr.stringcolumn(9).' ', \
                   SqR2arr=SqR2arr.stringcolumn(10).' ', \
                   EEJR2arr=EEJR2arr.stringcolumn(11).' ', $1):2

# width and height in inches; pixes per inch
w = 8.0
h = 6.0
ppi = 150

set term pngcairo enh col size (w*ppi),(h*ppi) font ",12"

nrow = 3
ncol = 1
l = 0.06
r = 0.06
t = 0.05
b = 0.06
dx = 0.08
dy = 0.07
load 'multiplot.cfg'
eval(init_margins(l, r, t, b, dx, dy, nrow, ncol))

stats fileSqLcorner;
nplot = STATS_records

load 'xyborder.cfg'
load 'grid.cfg'

do for [idx=0:nplot - 1] {

np = idx + 1

# Retrieve timestamp string
cmd = sprintf('time2str -t %s -z', word(tarr, np))
tstr = system(cmd)

# Retrieve longitude of equator crossing
phi = word(phiarr, np) + 0
phistr = sprintf('%.1f', phi)

# Retrieve LT of equator crossing
lt = word(ltarr, np) + 0
lthour = int(lt)
ltmin = int((lt - lthour) * 60)
ltstr = sprintf('%02d:%02d', lthour, ltmin)

# Retrieve satellite direction of equator crossing
dir = word(dirarr, np) + 0
if (dir > 0) {
  dirstr="upleg"
} else {
  dirstr="downleg"
}

cmd = sprintf('time2str -t %s -a -z', word(tarr, np))
fstr = system(cmd)
outstr = sprintf('%s/plot_%s.png', outdir, fstr)
set out outstr

str = sprintf('Generating plot %d/%d: %s...', np, nplot, outstr)
print str

set multiplot

set xrange [-60:60]
set xlabel "QD latitude (degrees)" offset 0,0.7
set ylabel "scalar residual (nT)" offset 1.3,0
load 'lines2.cfg'

eval(set_margins(1,1))

set title tstr.', {/Symbol \152} = '.phistr.'{/Symbol \260}'.', LT = '.ltstr.', kp = '.word(kparr, np).', '.dirstr.', track '.np

R2 = word(SqR2arr, np) + 0
R2str = sprintf('Sq R^2 = %.4f', R2)
set label 1 R2str at screen 0.82,0.72

plot fileF2 us 7:11 index idx w li lt 5 lw 4 ti "F^{(1)}", \
     fileF2 us 7:12 index idx w li lw 2 dt 2 ti "Sq internal", \
     fileF2 us 7:13 index idx w li lw 2 dt 2 ti "Sq external", \
     fileF2 us 7:($12 + $13) index idx w li lw 2 ti "Sq total", \
     fileF2 us 7:14 index idx w li lt 1 lw 4 ti "F^{(2)}"

unset label 1
unset title

set xrange [-30:30]

load 'xy2bordercol3.cfg'
set y2label "current density (mA/m)" offset -2,0

R2 = word(EEJR2arr, np) + 0
R2str = sprintf('EEJ R^2 = %.4f', R2)
set label 2 R2str at screen 0.82,0.4

eval(set_margins(2,1))
plot fileF2 us 7:14 index idx w li lw 4 ti "F^{(2)}", \
     fileF2 us 7:15 index idx w li lw 2 ti "F^{(2)} fit", \
     fileF2 us 7:(abs($7) <= 20 ? $16 : 1/0) index idx w li lw 4 ti "J_{/Symbol \146}" axes x1y2

unset label 2

nrow = 3
ncol = 2
eval(init_margins(l, r, t, b, dx, dy, nrow, ncol))

unset y2label
load 'xyborder.cfg'

set xlabel "residual norm ||y - A x||"
set ylabel "solution norm ||x||" offset 0,0
load 'lines.cfg'
load 'xylogon.cfg'

set xrange [1:1e4]

eval(set_margins(3,1))
plot fileSqLcurve us 2:3 index idx w lp lt 5 pt 7 ps 0.3 ti "Sq", \
     fileSqLcorner us 2:3 index idx w p lt 7 ps 3 pt 6 lw 4 ti ""

set xrange [*:*]
set yrange [*:*]
unset ylabel

eval(set_margins(3,2))
plot fileEEJLcurve us 2:3 index idx w lp lt 5 pt 7 ps 0.3 ti "EEJ", \
     fileEEJLcorner us 2:3 index idx w p lt 7 ps 3 pt 6 lw 4 ti ""

load 'xylogoff.cfg'
set yrange [*:*]

nrow = 3
ncol = 1
eval(init_margins(l, r, t, b, dx, dy, nrow, ncol))

unset multiplot

}
