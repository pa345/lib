#!/bin/sh
#
# Plot temporal histogram of data coverage for model
#
# NOTE: need to run mfield with -m option to generate datamap files

prefix_dir="output"
prefix="${prefix_dir}/datamap"

# Number of days for each histogram bin
binsize="30"

dt=$(echo "scale=10; ${binsize} / 365.25" | bc)

echo "prefix = ${prefix}"
echo "bin size = ${binsize} [days]"

tmin="2000"
tmax="2019"
nbins=$(echo "($tmax - $tmin) / $dt" | bc)

function proc
{
  file="$1"
  outfile="$2"

  #tmin=$(cat ${file} | datasel | awk '{print $1}' | head -1)
  #tmax=$(cat ${file} | datasel | awk '{print $1}' | tail -1)
  #nbins=$(echo "($tmax - $tmin) / $dt" | bc)

  echo "writing ${outfile}..."
  cat ${file} | datasel | awk '{print $1}' | gsl-histogram $tmin $tmax $nbins | awk '{print 0.5*($1+$2),$3}' > ${outfile}
}

rm -f ${prefix_dir}/hist_scal.sat*
rm -f ${prefix_dir}/hist_vec_Z.sat*

for sat in 0_1 $(seq 2 9); do
  fileF="${prefix}${sat}_F.dat"
  if [ -f ${fileF} ]; then
    echo "processing ${fileF}..."
    proc ${fileF} "${prefix_dir}/hist_scal.sat${sat}"
  fi

  #fileX="${prefix}${sat}_X.dat"
  #if [ -f ${fileX} ]; then
  #  echo "processing ${fileX}..."
  #  proc ${fileX} "hist_vec_X.sat${sat}"
  #fi

  #fileY="${prefix}${sat}_Y.dat"
  #if [ -f ${fileY} ]; then
  #  echo "processing ${fileY}..."
  #  proc ${fileY} "hist_vec_Y.sat${sat}"
  #fi

  fileZ="${prefix}${sat}_Z.dat"
  if [ -f ${fileZ} ]; then
    echo "processing ${fileZ}..."
    proc ${fileZ} "${prefix_dir}/hist_vec_Z.sat${sat}"
  fi

  fileZ_align="${prefix}${sat}_align_Z.dat"
  if [ -f ${fileZ_align} ]; then
    echo "processing ${fileZ_align}..."
    proc ${fileZ_align} "${prefix_dir}/hist_vec_Z_align.sat${sat}"
  fi

  #outfile_euler="hist_euler.sat${sat}"
  #echo "writing ${outfile_euler}..."
  #cat ${file} | datasel -c 10 --eq 1 | awk '{print $1}' | gsl-histogram $tmin $tmax $nbins > ${outfile_euler}
done

# now build final histogram file with all satellites suitable for gnuplot stacked histogram plotting
tmpfile=$(mktemp)

outfile="${prefix_dir}/hist_vec_Z.txt"
paste ${prefix_dir}/hist_vec_Z.sat* > ${tmpfile}
echo "writing $outfile"
cat $tmpfile | awk '{print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18}' > $outfile

outfile="${prefix_dir}/hist_vec_Z_align.txt"
paste ${prefix_dir}/hist_vec_Z_align.sat* > ${tmpfile}
echo "writing $outfile"
cat $tmpfile | awk '{print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18}' > $outfile

outfile="${prefix_dir}/hist_scal.txt"
paste ${prefix_dir}/hist_scal.sat* > ${tmpfile}
echo "writing $outfile"
cat $tmpfile | awk '{print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18}' > $outfile

rm -f $tmpfile
