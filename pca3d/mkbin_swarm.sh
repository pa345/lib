#!/bin/sh
#
# Read all Sswarm CDF files and process into one magdata formatted file

cfgfile="INVERT_preproc.cfg"

# Generate latest Swarm index files
#sh ./mkidx_swarm.sh

current_year=$(date +"%Y")

idxfile=$(mktemp)

for sat in A B; do
  satfile="./swarm${sat}.idx"
  outfile="data_invert/swarm${sat}.dat"
  rm -f ${outfile}

  for year in $(seq 2013 ${current_year}); do

    grep ".*MAG${sat}_LR_1B_${year}.*.cdf" ${satfile} > ${idxfile}

    echo "Processing satellite ${sat} for year ${year}..."
    if [ "${year}" = "2013" ]; then
      ./invert_preproc -s ${idxfile} -C ${cfgfile} -N "Swarm ${sat}" -o ${outfile}
    else
      ./invert_preproc -s ${idxfile} -C ${cfgfile} -N "Swarm ${sat}" --append ${outfile} -o ${outfile}
    fi
  done
done

rm -f $idxfile
