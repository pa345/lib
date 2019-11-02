#!/bin/sh
#
# Read all Cryosat CDF files and process into one magdata formatted file

cfgfile="MF_cryosat_preproc.cfg"

for satnum in $(seq 1 3); do
  datadir="$DATAHOME/Cryosat/Stage1_FGM${satnum}"
  outfile="data/cryosat${satnum}.dat"

  rm -f ${outfile}

  for year in $(seq 2010 2018); do
    idxfile=$(mktemp)
    echo "Generating index file for FGM${satnum} for year ${year}..."
    find $datadir/$year -name "*.cdf" | sort -g > ${idxfile}

    echo "Processing year ${year}..."
    if [ "${year}" = "2010" ]; then
      ./mfield_preproc -r ${idxfile} -C ${cfgfile} -o ${outfile}
    else
      ./mfield_preproc -r ${idxfile} -C ${cfgfile} --append ${outfile} -o ${outfile}
    fi

    rm -f $idxfile
  done
done