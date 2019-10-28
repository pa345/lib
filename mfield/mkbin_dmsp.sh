#!/bin/sh
#
# Read all DMSP CDF files and process into one magdata formatted file

datadir="$DATAHOME/DMSP/MAG/Stage2"
cfgfile="MF_dmsp_preproc.cfg"

for sat in f16 f17 f18; do
  outfile="data/${sat}.dat"
  rm -f ${outfile}

  for year in $(seq 2009 2016); do
    datafile="${datadir}/${sat}_${year}_Stage2.cdf"

    echo "Processing year ${year}..."
    if [ "${year}" = "2009" ]; then
      ./mfield_preproc -D ${datafile} -C ${cfgfile} -o ${outfile}
    else
      ./mfield_preproc -D ${datafile} -C ${cfgfile} --append ${outfile} -o ${outfile}
    fi
  done
done
