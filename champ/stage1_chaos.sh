#!/bin/bash
#
# This script computes main,crustal,external field values for CHAMP data
#
# Usage: ./stage1.sh [year]

# load environment variables
. $HOME/.bashrc

cdfdir="/data/CHAMP/L3"
outdir="/data/CHAMP/Stage1_CHAOS"

prog="/home/palken/usr/bin/stage1"

# Use CHAOS external field model?
extra_flags="-h"

year=""
if test -n "$1"; then
  year="$1"
fi

#rm -rf ${outdir}
mkdir ${outdir}

if [ -z "${year}" ]; then
  yearlist=$(seq 2000 2010)
else
  yearlist="$year"
fi

for yr in ${yearlist}; do
  filelist=$(ls ${cdfdir}/${yr}/*.cdf)
#  filelist=$(ls ${cdfdir}/*.cdf)

  for file in ${filelist}; do
    bname=$(basename "${file}" .cdf)
    outfile="${outdir}/${bname}_Stage1.cdf"

    # Check first if we already processed this file
    if [[ ! -f "${outfile}" ]]; then
      echo "Processing: ${file}"
      ${prog} -c ${file} -o ${outfile} ${extra_flags}
    fi
  done
done
