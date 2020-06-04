#!/bin/sh
#
# Grid the original data, model, and residuals for easy plotting

output_dir="output_Swarm"
prefix="swarmA"
satnum="0"

outfile="${prefix}_mod_X.grd"
echo "Generating ${outfile}..."
cat ${output_dir}/res${satnum}_X_iter1.dat | gridgeo -a 0 -b 24 -c -90 -d 90 -x 2 -y 5 -z 12 -u 1 -v 1 -o ${outfile}

outfile="${prefix}_mod_Y.grd"
echo "Generating ${outfile}..."
cat ${output_dir}/res${satnum}_Y_iter1.dat | gridgeo -a 0 -b 24 -c -90 -d 90 -x 2 -y 5 -z 12 -u 1 -v 1 -o ${outfile}

outfile="${prefix}_mod_Z.grd"
echo "Generating ${outfile}..."
cat ${output_dir}/res${satnum}_Z_iter1.dat | gridgeo -a 0 -b 24 -c -90 -d 90 -x 2 -y 5 -z 12 -u 1 -v 1 -o ${outfile}

outfile="${prefix}_res_X.grd"
echo "Generating ${outfile}..."
cat ${output_dir}/res${satnum}_X_iter1.dat | gridgeo -a 0 -b 24 -c -90 -d 90 -x 2 -y 5 -z 13 -u 1 -v 1 -o ${outfile}

outfile="${prefix}_res_Y.grd"
echo "Generating ${outfile}..."
cat ${output_dir}/res${satnum}_Y_iter1.dat | gridgeo -a 0 -b 24 -c -90 -d 90 -x 2 -y 5 -z 13 -u 1 -v 1 -o ${outfile}

outfile="${prefix}_res_Z.grd"
echo "Generating ${outfile}..."
cat ${output_dir}/res${satnum}_Z_iter1.dat | gridgeo -a 0 -b 24 -c -90 -d 90 -x 2 -y 5 -z 13 -u 1 -v 1 -o ${outfile}

outfile="${prefix}_data_X.grd"
echo "Generating ${outfile}..."
cat ${output_dir}/res${satnum}_X_iter1.dat | datasel | awk '{print $2,$5,$10-$11}' | gridgeo -a 0 -b 24 -c -90 -d 90 -u 1 -v 1 -o ${outfile}

outfile="${prefix}_data_Y.grd"
echo "Generating ${outfile}..."
cat ${output_dir}/res${satnum}_Y_iter1.dat | datasel | awk '{print $2,$5,$10-$11}' | gridgeo -a 0 -b 24 -c -90 -d 90 -u 1 -v 1 -o ${outfile}

outfile="${prefix}_data_Z.grd"
echo "Generating ${outfile}..."
cat ${output_dir}/res${satnum}_Z_iter1.dat | datasel | awk '{print $2,$5,$10-$11}' | gridgeo -a 0 -b 24 -c -90 -d 90 -u 1 -v 1 -o ${outfile}

