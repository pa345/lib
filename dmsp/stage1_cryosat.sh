#!/bin/sh
#
# This script computes main field values for Cryosat data
#
# Usage: stage1.sh [year]

# load environment variables
#. $HOME/.bashrc
export LD_LIBRARY_PATH="/home/palken/usr/lib:/usr/local/MATLAB/MATLAB_Runtime/v91/runtime/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v91/bin/glnxa64:/usr/local/MATLAB/MATLAB_Runtime/v91/sys/os/glnxa64:/usr/lib64:/usr/local/cdf/lib:"

# year to process
year="2009"
if test -n "$1"; then
  year="$1"
fi

indir="$DATAHOME/Cryosat/Original_Data"
outdir="$DATAHOME/Cryosat/Stage1_FGM3"

prog="$HOME/usr/bin/stage1"

# Extra flags (such as use CHAOS)
# Aug 6 2018: CHAOS-6 has a bug in the external field coefficients for years 2011-2013 due
# to lack of satellite data
#extra_flags="-h"

echo "==== PROCESSING YEAR ${year} ===="

mkdir -p $outdir/$year

for file in $(ls ${indir}/CS_OPER_*_${year}*.cdf); do
  bname=$(basename "${outdir}/${file}" _Stage0.cdf)
  outfile="${outdir}/${year}/${bname}_Stage1.cdf"

  # Check first if we already processed this file
  if [[ ! -f "${outfile}" ]]; then
    echo "Processing: ${file}"
    ${prog} -r ${file} -o ${outfile} ${extra_flags}
  fi
done
