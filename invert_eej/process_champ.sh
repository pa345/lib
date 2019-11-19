#!/bin/sh
#
# Run EEF chain on CHAMP data

prog="./main"

# set for profiles only
flags="-p -a 0.00001 -b 23.99999"
#flags=""

for year in $(seq -w 10 10); do
  log_dir="log_champ_all_${year}"
  idx_file="champ${year}.idx"

  rm -rf ${log_dir}
  mkdir ${log_dir}

  echo "Processing year ${year} (${idx_file}, ${log_dir})"
  screen -d -m -S champ${year} ${prog} -c ${idx_file} -l ${log_dir} ${flags}
done
