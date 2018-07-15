#!/bin/sh
#
# Make index files for Swarm, using the latest baseline
# version available for each day

datadir="$DATAHOME/SWARM/MAG/Stage1"

function dosat
{
  sat_name="$1"
  idxfile="$2"

  rm -f ${idxfile}

  echo "sat = ${sat_name}"

  # loop over years
  for year in $(seq 2013 2018); do
    echo "year = ${year}"

    # loop over months
    for month in $(seq -w 01 12); do

      # loop over days
      for day in $(seq -w 01 31); do

        # check for duplicate files for the same date, and choose
        # the latest one (ie with the highest version number);
        # also require a minimum baseline of 03
        vermax="0300"
        filemax=""
        for file in $(find ${datadir} -name "SW*MAG${sat_name}_LR_1B_${year}${month}${day}*.cdf"); do
          version=$(basename ${file} | awk -F'_' '{print $8}')
          if [ "${version}" -gt "${vermax}" ]; then
            vermax=${version}
            filemax=${file}
          fi
        done

        if [ -n "${filemax}" ]; then
          echo "${filemax}" >> ${idxfile}
        fi
      done
    done
  done

  echo "Output is ${idxfile}"
}

dosat "A" "swarmA.idx"
dosat "B" "swarmB.idx"
dosat "C" "swarmC.idx"
