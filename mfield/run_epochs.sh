#!/bin/sh
#
# Run magnetic field model on all epochs

start_epoch="2009.0"
end_epoch="2016.0"
numit="5"

input_files="data/champ*.dat data/F17a.dat data/F17.dat data/swarm*.dat data/obs/*SV.dat"

prog="./mfield"
coef_dir="$MYLIBHOME/mfield/coef_GAP"
cfgfile="MF_dmsp.cfg"

# SV/SA damping parameter
lambda_sv="0.095"
lambda_sa="1.0"

rm -rf ${coef_dir}
mkdir ${coef_dir}

# step size = 1 month = 1/12 years
step="0.083333"

# data window size in years
window_size="3"

#extra_flags="-r -d -p 30"
extra_flags="-p 30"

for epoch in $(seq ${start_epoch} ${step} ${end_epoch}); do

  # find 3-year window centered on epoch
  low=$(echo ${epoch} - 0.5*${window_size} | bc)
  high=$(echo ${epoch} + 0.5*${window_size} | bc)

  echo "epoch = ${epoch} low = ${low} high = ${high} window_size = ${window_size}"

  coef_file="${coef_dir}/coef.${epoch}.txt"
  echo "computing model for epoch ${epoch}; output in ${coef_file}..."
  ${prog} -n ${numit} -e ${epoch} -C ${cfgfile} -v ${lambda_sv} -a ${lambda_sa} --tmin ${low} --tmax ${high} -o ${coef_file} ${extra_flags} ${input_files}
done
