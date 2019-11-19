#!/bin/sh
#
# Make a GMT scatter plot of (x,y) data
# Usage: gmtscat.sh [xyfile]

infile="datamap.dat"
if test -n "$1"; then
  infile="$1"
fi

outfile="${infile}.ps"
rm -f $outfile

echo "input file = ${infile}"

xyfile=$(mktemp)
cat ${infile} | datasel | awk '{print $2,$3}' > ${xyfile}

# size of each scattered data point (circle)
pointsize="0.005i"

# Equatorial view
gmt pscoast -X0.8i -Y6.5i -R-180/180/-90/90 -JW0/6.8 -B45g0f0/10g0f0 -N1 -Di -W2 -A10/1 -G180/180/180 -V -P -K >> $outfile
gmt psxy $xyfile -: -R -JW -O -P -G255/0/0 -Sc${pointsize} -V -K >> $outfile

# North pole
gmt pscoast -X0.4i -Y-2.6i -R-180/180/40/90 -JG0/90/2 -B45g0f0/10g0f0 -N1 -Di -W2 -A10/1 -G180/180/180 -V -P -O -K >> $outfile
gmt psxy $xyfile -: -R -JG -O -P -G255/0/0 -Sc${pointsize} -V -K >> $outfile

# South pole
gmt pscoast -X4.0i -Y0 -R-180/180/-90/-40 -JG0/-90/2 -B45g0f0/10g0f0 -N1 -Di -W2 -A10/1 -G180/180/180 -V -P -O -K >> $outfile
gmt psxy $xyfile -: -R -JG -O -P -G255/0/0 -Sc${pointsize} -V >> $outfile

rm -f ${xyfile}

echo "Output is ${outfile}"
