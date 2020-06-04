#!/bin/sh
#
# Read all Cryosat CDF files and process into one magdata formatted file

datadir="$DATAHOME/Cryosat/Stage1_FGM1"

outfile="data_invert/cryosat.dat"
cfgfile="INVERT_preproc.cfg"

rm -f ${outfile}

idxfile=$(mktemp)

for year in $(seq 2010 2018); do
  echo "Generating index file for year ${year}..."
  find $datadir -name "CS_OPER_MAG_${year}*.cdf" | sort -g > ${idxfile}

  echo "Processing year ${year}..."
  if [ "${year}" = "2010" ]; then
    ./invert_preproc -r ${idxfile} -C ${cfgfile} -N "Cryosat FGM1" -o ${outfile}
  else
    ./invert_preproc -r ${idxfile} -C ${cfgfile} -N "Cryosat FGM1" --append ${outfile} -o ${outfile}
  fi
done

rm -f $idxfile
