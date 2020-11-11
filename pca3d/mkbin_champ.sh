#!/bin/sh
#
# Read all CHAMP CDF files and process into one magdata formatted file

datadir="$DATAHOME/CHAMP/Stage1_CHAOS"
rmsfile="CHAMP_rms.txt"

outfile="data_invert/champ.dat"
cfgfile="INVERT_preproc.cfg"

rm -f ${outfile} ${rmsfile}

idxfile=$(mktemp)
tmpfile=$(mktemp)

for year in $(seq 2000 2010); do
  echo "Generating index file for year ${year}..."
  find $datadir -name "CH*${year}*.cdf" | sort -g > ${idxfile}

  echo "Processing year ${year}..."
  if [ "${year}" = "2000" ]; then
    ./invert_preproc -c ${idxfile} -C ${cfgfile} -N "CHAMP" --rms_file ${tmpfile} -o ${outfile}
  else
    ./invert_preproc -c ${idxfile} -C ${cfgfile} -N "CHAMP" --rms_file ${tmpfile} --append ${outfile} -o ${outfile}
  fi

  cat ${tmpfile} >> ${rmsfile}
done

rm -f $idxfile $tmpfile
