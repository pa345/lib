#!/bin/sh

coefdir1="maps_CHAOS"
coefdir2="maps_G"
coefdir3="maps_D"
coefdir="maps_combined"

mkdir ${coefdir}

idx=1
for epoch in $(seq 2000.0 0.1 2019.0); do
  infile1="${coefdir1}/map.${epoch}.png"
  infile2="${coefdir2}/map.${epoch}.png"
  infile3="${coefdir3}/map.${epoch}.png"

  istr=$(seq -f "%03g" $idx $idx)
  outfile="${coefdir}/map.${istr}.png"

  echo "generating map for epoch ${epoch} (${outfile})..."
  #convert ${infile1} ${infile2} ${infile3} +append ${outfile}
  convert ${infile1} ${infile2} +append ${outfile}

  idx=$((idx+1))
done

#ffmpeg -framerate 4 -i map.%03d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -vf fps=25 ${outfile}
