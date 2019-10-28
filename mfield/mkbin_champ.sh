#!/bin/sh
#
# Read all CHAMP CDF files and process into one magdata formatted file

datadir="$DATAHOME/CHAMP/Stage1_CHAOS"
outfile="data/champ.dat"

rm -f ${outfile}

for year in $(seq 2000 2010); do
  idxfile=$(mktemp)
  echo "Generating index file for year ${year}..."
  find $datadir -name "CH*${year}*.cdf" | sort -g > ${idxfile}

  echo "Processing year ${year}..."
  if [ "${year}" = "2000" ]; then
    ./mfield_preproc -c ${idxfile} -C MF_preproc.cfg -o ${outfile}
  else
    ./mfield_preproc -c ${idxfile} -C MF_preproc.cfg --append ${outfile} -o ${outfile}
  fi

  rm -f $idxfile
done
