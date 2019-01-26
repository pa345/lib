#!/bin/sh
#
# This script computes main field values for DMSP data
#
# Usage: stage1.sh [FXX] [year]

# load environment variables
. $HOME/.bashrc

# satellite to process
dmsp_sat="F15"
if test -n "$1"; then
  dmsp_sat="$1"
fi

# year to process
year="2009"
if test -n "$2"; then
  year="$2"
fi

indir="$DATAHOME/DMSP/MAG/Stage0"
outdir="$DATAHOME/DMSP/MAG/Stage1"

prog="$HOME/usr/bin/stage1"

# Extra flags (such as use CHAOS)
# Aug 6 2018: CHAOS-6 has a bug in the external field coefficients for years 2011-2013 due
# to lack of satellite data
#extra_flags="-h"

echo "==== PROCESSING YEAR ${year} ===="

mkdir -p $outdir/$year

for file in $(ls ${indir}/${year}/DMSP_${dmsp_sat}_*.cdf); do
  bname=$(basename "${outdir}/${year}/${file}" _Stage0.cdf)
  outfile="${outdir}/${year}/${bname}_Stage1.cdf"

  # Check first if we already processed this file
  if [[ ! -f "${outfile}" ]]; then
    echo "Processing: ${file}"
    ${prog} -d ${file} -o ${outfile} ${extra_flags}
  fi
done
