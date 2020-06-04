#!/bin/sh
#
# Read all CHAMP CDF files and process into one magdata formatted file

datadir="$DATAHOME/CHAMP/Stage1_CHAOS"

outfile="data_invert/champ.dat"
cfgfile="INVERT_preproc.cfg"

rm -f ${outfile}

idxfile=$(mktemp)

for year in $(seq 2000 2010); do
  echo "Generating index file for year ${year}..."
  find $datadir -name "CH*${year}*.cdf" | sort -g > ${idxfile}

  echo "Processing year ${year}..."
  if [ "${year}" = "2000" ]; then
    ./invert_preproc -c ${idxfile} -C ${cfgfile} -N "CHAMP" -o ${outfile}
  else
    ./invert_preproc -c ${idxfile} -C ${cfgfile} -N "CHAMP" --append ${outfile} -o ${outfile}
  fi
done

rm -f $idxfile
