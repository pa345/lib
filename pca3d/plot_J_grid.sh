#!/bin/sh

data_dir="./output_350km"

for file in $(ls ${data_dir}/J_grid_*.txt); do
  ts=$(head ${file} | grep Timestamp | sed -e 's/# Timestamp: \([0-9]*\)/\1/')
  tstr=$(time2str -t ${ts})
  idxnum=$(echo ${file} | sed -e 's/.*J_grid_\(.*\).txt/\1/')
  #~/usr/bin/gnuplot -e "idxnum='${idxnum}'; tstr='${tstr}'" ./plot_J_grid.gp

  bname=$(basename ${file} .txt)
  outfile="${data_dir}/${bname}.png"
  title="${tstr}"

  # 110 km J data
  # J_r
  #/data/palken/mainlib/msynth/plots/genmap.sh -i ${file} -o ${outfile} --cbmin -0.1 --cbmax 0.1 --cbstep 0.05 --cblabel "uA/m^2" -z 3 -t "${title}" --xnodes 181 --ynodes 91
  # J_t
  #/data/palken/mainlib/msynth/plots/genmap.sh -i ${file} -o ${outfile} --cbmin -1.5 --cbmax 1.5 --cbstep 0.5 --cblabel "uA/m^2" -z 4 -t "${title}" --xnodes 181 --ynodes 91

  # 350 km J data
  tmpfile=$(mktemp)
  cat ${file} | datasel | awk '{print $1,$2,$5*1e3}' > ${tmpfile}
  /data/palken/mainlib/msynth/plots/genmap.sh -i ${tmpfile} -o ${outfile} --cbmin -1 --cbmax 1 --cbstep 0.5 --cblabel "nA/m^2" -z 3 -t "${title}" --xnodes 181 --ynodes 91
  #rm -f ${tmpfile}
done

rm -f ${data_dir}/*.eps ${data_dir}/*.ps

#ffmpeg -framerate 4 -pattern_type glob -i '*.png' -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -vf fps=25 ${outfile}
