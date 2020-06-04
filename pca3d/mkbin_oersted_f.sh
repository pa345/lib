#!/bin/sh
#
# Read all Oersted (f) CDF files and process into one magdata formatted file

datadir="$DATAHOME/Oersted/MAG/Stage1_scalar"

outfile="data_invert/oersted_scalar.dat"
cfgfile="INVERT_preproc.cfg"

rm -f ${outfile}

idxfile=$(mktemp)

for year in $(seq 1999 2013); do
  echo "Generating index file for year ${year}..."
  find $datadir -name "Oersted_${year}*.cdf" | sort -g > ${idxfile}

  echo "Processing year ${year}..."
  if [ "${year}" = "1999" ]; then
    ./invert_preproc -F ${idxfile} -C ${cfgfile} -N "Oersted scalar" -o ${outfile}
  else
    ./invert_preproc -F ${idxfile} -C ${cfgfile} -N "Oersted scalar" --append ${outfile} -o ${outfile}
  fi
done

rm -f $idxfile
