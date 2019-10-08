#!/bin/sh
#
# This script runs stage2a on Stage1 DMSP data files

dmsp_sats="f16 f17 f18"
prog="./stage2a"
out_dir="/data/DMSP/MAG/Stage1a"

mkdir -p ${out_dir}

for sat in ${dmsp_sats}; do
  echo $sat
  for year in $(seq 2009 2016); do
    idxfile="./${sat}_${year}.idx"
    outfile="${out_dir}/${sat}_${year}_Stage1a.cdf"
    #eulerfile="./euler_${sat}.txt"
    #${prog} -i ${idxfile} -e ${eulerfile} -o ${outfile}
    ${prog} -i ${idxfile} -o ${outfile}
  done
done
