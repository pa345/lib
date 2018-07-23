#!/bin/sh
#
# Make index file for CHAMP

datadir="$DATAHOME/CHAMP/Stage1_CHAOS"

for year in $(seq 2000 2010); do
  idxfile="champ_${year}.idx"
  echo "Generating $idxfile..."
  find $datadir -name "CH*${year}*.cdf" | sort -g > ${idxfile}
done
