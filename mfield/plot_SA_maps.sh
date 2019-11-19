#!/bin/sh

set term pngcairo enh col size 1000,1000

mapprog="$DATAHOME/palken/repo/msynth/src/print_map"
plotprog="$DATAHOME/palken/repo/msynth/src/plots/genmap.sh"

#coef_file="/data/palken/repo/msynth/src/cof/CHAOS-6-x8_core.shc"
#coefdir="maps_CHAOS"
#title="CHAOS-6-x8"

coef_file="Model_AAA3.shc"
coefdir="maps_AAA3"
title="Model A"

plot_args="-c "uT/yr^2" --cbmin -2.0 --cbmax 2.0 --cbstep 1.0"

start_time="2000.6"
end_time="2020.0"
time_step="0.1"

# maximum SH degree for SA maps
nmax="8"

mkdir ${coefdir}

idx=1
for epoch in $(seq ${start_time} ${time_step} ${end_time}); do
  outfile="${coefdir}/map.${epoch}.png"

  echo "generating SA map for epoch ${epoch}..."
  tmpfile=$(mktemp)
  ${mapprog} -h ${coef_file} -e ${epoch} -n $nmax -o $tmpfile
  ${plotprog} -i $tmpfile -o ${outfile} -t "${title}, epoch ${epoch}" ${plot_args}
  rm -f $tmpfile
done

# cleanup
rm -f ${coefdir}/*.eps ${coefdir}/*.ps

#ffmpeg -framerate 4 -pattern_type glob -i '*.png' -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -vf fps=25 ${outfile}
