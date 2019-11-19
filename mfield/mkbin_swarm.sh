#!/bin/sh
#
# Read all Sswarm CDF files and process into one magdata formatted file

cfgfile="MF_preproc.cfg"
#suffix="_nocrust"
suffix=""

# Generate latest Swarm index files
#sh ./mkidx_swarm.sh

current_year=$(date +"%Y")

for sat in A B; do
  satfile="./swarm${sat}.idx"
  outfile="data/swarm${sat}${suffix}.dat"
  rm -f ${outfile}

  for year in $(seq 2013 ${current_year}); do
    idxfile=$(mktemp)

    grep ".*MAG${sat}_LR_1B_${year}.*.cdf" ${satfile} > ${idxfile}

    echo "Processing satellite ${sat} for year ${year}..."
    if [ "${year}" = "2013" ]; then
      ./mfield_preproc -s ${idxfile} -C ${cfgfile} -N "Swarm ${sat}" -o ${outfile}
    else
      ./mfield_preproc -s ${idxfile} -C ${cfgfile} -N "Swarm ${sat}" --append ${outfile} -o ${outfile}
    fi

    rm -f $idxfile
  done
done
