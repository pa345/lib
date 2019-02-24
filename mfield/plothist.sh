#!/bin/sh
#
# Plot temporal histogram of data coverage for model
#
# NOTE: need to run mfield with -m option to generate datamap files

prefix="output/datamap"

# Number of days for each histogram bin
binsize="10"

dt=$(echo "scale=10; ${binsize} / 365.25" | bc)

echo "prefix = ${prefix}"
echo "bin size = ${binsize} [days]"

function proc
{
  file="$1"
  outfile="$2"

  tmin=$(cat ${file} | datasel | awk '{print $1}' | head -1)
  tmax=$(cat ${file} | datasel | awk '{print $1}' | tail -1)

  nbins=$(echo "($tmax - $tmin) / $dt" | bc)

  echo "writing ${outfile}..."
  cat ${file} | datasel | awk '{print $1}' | gsl-histogram $tmin $tmax $nbins > ${outfile}
}

for sat in $(seq 0 1); do
  fileF="${prefix}${sat}_F.dat"
  if [ -f ${fileF} ]; then
    echo "processing ${fileF}..."
    proc ${fileF} "hist_scal.sat${sat}"
  fi

  fileX="${prefix}${sat}_X.dat"
  if [ -f ${fileX} ]; then
    echo "processing ${fileX}..."
    proc ${fileX} "hist_vec_X.sat${sat}"
  fi

  fileY="${prefix}${sat}_Y.dat"
  if [ -f ${fileY} ]; then
    echo "processing ${fileY}..."
    proc ${fileY} "hist_vec_Y.sat${sat}"
  fi

  fileZ="${prefix}${sat}_Z.dat"
  if [ -f ${fileZ} ]; then
    echo "processing ${fileZ}..."
    proc ${fileZ} "hist_vec_Z.sat${sat}"
  fi

  #outfile_euler="hist_euler.sat${sat}"
  #echo "writing ${outfile_euler}..."
  #cat ${file} | datasel -c 10 --eq 1 | awk '{print $1}' | gsl-histogram $tmin $tmax $nbins > ${outfile_euler}
done
