#!/home/palken/usr/bin/gnuplot
#
# Plot F2 curves and corresponding L-curves for Sq fit
#
# Adjust 'logdir' and 'outdir' below

logdir = 'log'
outdir = 'log/plots'

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
plot fileprof us ( tarr=tarr.stringcolumn(2).' ', \
                   phiarr=phiarr.stringcolumn(3).' ', \
                   ltarr=ltarr.stringcolumn(4).' ', \
                   kparr=kparr.stringcolumn(8).' ', \
                   dirarr=dirarr.stringcolumn(9).' ', $1):2

nrow = 3
ncol = 1

load 'multi_default.cfg'

l = 0.5
r = 0.2
b = 0.1
vbuffer = 0.5
plotwidth = 5.0
plotheight = 1.2
fontsize = ",12"

load 'multi_defs.cfg'
load 'multi_png.cfg'

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

set multiplot layout nrow,ncol

load 'multi_reset.cfg'

set xrange [-60:60]
set xlabel "QD latitude (degrees)"
set ylabel "scalar residual (nT)"
load 'lines2.cfg'

set title tstr.', {/Symbol \152} = '.phistr.'{/Symbol \260}'.', LT = '.ltstr.', kp = '.word(kparr, np).', '.dirstr.', track '.np

plot fileF2 us 7:11 index idx w li lt 5 lw 4 ti "F^{(1)}", \
     fileF2 us 7:12 index idx w li lw 2 dt 2 ti "Sq internal", \
     fileF2 us 7:13 index idx w li lw 2 dt 2 ti "Sq external", \
     fileF2 us 7:($12 + $13) index idx w li lw 2 ti "Sq total", \
     fileF2 us 7:14 index idx w li lt 1 lw 4 ti "F^{(2)}"

unset title

load 'incrow.cfg'

set xrange [-30:30]
#set y2range [-0.25:0.15]

load 'xy2bordercol3.cfg'
set y2label "current density (mA/m)"

plot fileF2 us 7:14 index idx w li lw 4 ti "F^{(2)}", \
     fileF2 us 7:15 index idx w li lw 2 ti "F^{(2)} fit", \
     fileEEJ us 1:2 index idx w li lw 4 ti "J_{/Symbol \146}" axes x1y2

unset y2label
load 'xyborder.cfg'

load 'incrow.cfg'

set xlabel "residual norm ||y - A x||"
set ylabel "solution norm ||x||"
load 'lines.cfg'
load 'xylogon.cfg'

set xrange [1:1e4]
#set yrange [.1:1e6]

set rmargin at screen 0.45

plot fileSqLcurve us 2:3 index idx w lp lt 5 pt 7 ps 0.3 ti "Sq", \
     fileSqLcorner us 2:3 index idx w p lt 7 ps 3 pt 6 lw 4 ti ""

set lmargin at screen 0.55
set rmargin at screen right(currentcolumn,ncol,w,l,r,hbuffer)
#set xrange [1e-2:1e3]
set xrange [*:*]
#set yrange [1e-4:1e10]
set yrange [*:*]
unset ylabel

plot fileEEJLcurve us 2:3 index idx w lp lt 5 pt 7 ps 0.3 ti "EEJ", \
     fileEEJLcorner us 2:3 index idx w p lt 7 ps 3 pt 6 lw 4 ti ""

load 'xylogoff.cfg'
set yrange [*:*]

unset multiplot

}
