#!/bin/sh
#
# This script computes main, crustal, and external field
# values for Swarm data
#
# It is designed to be run continually, so as new data
# files are downloaded, the script checks if a file has
# already been processed before computing the fields
#
# Usage: ./stage1.sh [A|B|C] [year]

export LD_LIBRARY_PATH="/home/palken/usr/lib:/usr/local/MATLAB/MATLAB_Runtime/v91/runtime/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v91/bin/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v91/sys/os/glnxa64:/usr/lib64:/usr/local/cdf/lib:"

cdfdir="$DATAHOME/SWARM/MAG/Unzipped_Data"
outdir="$DATAHOME/SWARM/MAG/Stage1"
lpdir="$DATAHOME/SWARM/EFI/LP_Unzipped"

# Use CHAOS external field model
extra_flags="-h"
#extra_flags=""

prog="/home/palken/usr/bin/stage1"

sat=""
if test -n "$1"; then
  sat="$1"
fi

year=""
if test -n "$2"; then
  year="$2"
fi

#rm -rf ${outdir}
mkdir ${outdir}

for file in $(ls ${cdfdir}/*MAG${sat}*LR_1B_${year}*.cdf); do
  bname=$(basename "${file}" .cdf)
  outfile="${outdir}/${bname}_Stage1.cdf"

  # Find LP data file for this date
  lp_arg=""
  cursat=$(echo ${bname} | sed -r 's/.*MAG([A-Z]).*/\1/')
  tstr=$(echo ${bname} | sed -r 's/.*_1B_([0-9]*).*/\1/')
  lpfile=$(find ${lpdir} -name "SW*EFI${cursat}*_1B_${tstr}T*.txt")
  if [[ -f "${lpfile}" ]]; then
    lp_arg="-l ${lpfile}"
  fi

  # Check first if we already processed this file
  if [[ ! -f "${outfile}" ]]; then
    echo "Processing: ${file}"
    ${prog} -s ${file} -o ${outfile} ${lp_arg} ${extra_flags}
  fi
done
