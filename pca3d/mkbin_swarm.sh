#!/bin/sh
#
# Read all Swarm CDF files and process into one magdata formatted file

cfgfile="INVERT_preproc.cfg"

# Generate latest Swarm index files
#sh ./mkidx_swarm.sh

current_year=$(date +"%Y")

idxfile=$(mktemp)
tmpfile=$(mktemp)

for sat in A B; do
  satfile="./swarm${sat}.idx"
  outfile="data_invert/swarm${sat}.dat"
  rmsfile="swarm${sat}_rms.txt"

  rm -f ${outfile} ${rmsfile}

  for year in $(seq 2013 ${current_year}); do

    grep ".*MAG${sat}_LR_1B_${year}.*.cdf" ${satfile} > ${idxfile}

    echo "Processing satellite ${sat} for year ${year}..."
    if [ "${year}" = "2013" ]; then
      ./invert_preproc -s ${idxfile} -C ${cfgfile} -N "Swarm ${sat}" --rms_file ${tmpfile} -o ${outfile}
    else
      ./invert_preproc -s ${idxfile} -C ${cfgfile} -N "Swarm ${sat}" --rms_file ${tmpfile} --append ${outfile} -o ${outfile}
    fi

  cat ${tmpfile} >> ${rmsfile}
  done
done

rm -f $idxfile $tmpfile
