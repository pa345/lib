#!/bin/sh
#
# Read all Oersted CDF files and process into one magdata formatted file

datadir="$DATAHOME/Oersted/MAG/Stage1_vector"

outfile="data_invert/oersted_vector.dat"
cfgfile="INVERT_preproc.cfg"

rm -f ${outfile}

idxfile=$(mktemp)

for year in $(seq 1999 2005); do
  echo "Generating index file for year ${year}..."
  find $datadir -name "Oersted_${year}*.cdf" | sort -g > ${idxfile}

  echo "Processing year ${year}..."
  if [ "${year}" = "1999" ]; then
    ./invert_preproc -L ${idxfile} -C ${cfgfile} -N "Oersted vector" -o ${outfile}
  else
    ./invert_preproc -L ${idxfile} -C ${cfgfile} -N "Oersted vector" --append ${outfile} -o ${outfile}
  fi
done

rm -f $idxfile
