#!/bin/sh

band="11"
alt="110"

cbmin="-0.3"
cbmax="0.3"

data_dir="modes/${alt}/band_${band}"

for mode in $(seq -w 1 10); do
  echo "Generating ${alt} km mode ${mode}..."
  sh /data/palken/mainlib/msynth/plots/genmap.sh -i ${data_dir}/J_mode_${mode}.txt -z 5 --cbmin ${cbmin} --cbmax ${cbmax} --no-colorbar --xnodes 360 --ynodes 180 -t "J@-@~\152@~@- real mode ${mode}, ${alt} km" -o ${data_dir}/J_${mode}_phi_re.png
  sh /data/palken/mainlib/msynth/plots/genmap.sh -i ${data_dir}/J_mode_${mode}.txt -z 8 --cbmin ${cbmin} --cbmax ${cbmax} --no-colorbar --xnodes 360 --ynodes 180 -t "J@-@~\152@~@- imag mode ${mode}, ${alt} km" -o ${data_dir}/J_${mode}_phi_im.png
done

exit

cd modes/110/1_cpd
sh /data/palken/mainlib/msynth/plots/genmap.sh -i J_mode_01.txt -z 8 --cbmin -600 --cbmax 600 --cbstep 300 --xnodes 360 --ynodes 180 --no-colorbar -t "J@-@~\152@~@- spatial mode 1, 110km" -o J_01_phi_110km.png
sh /data/palken/mainlib/msynth/plots/genmap.sh -i J_mode_05.txt -z 3 --cbmin -50 --cbmax 50 --cbstep 25 --xnodes 360 --ynodes 180 --no-colorbar -t "J@-r@- spatial mode 5, 110km" -o J_05_r_110km.png
convert +append J_01_phi_110km.png J_05_r_110km.png out_110.png

#cd modes/300/1_cpd
#sh /data/palken/mainlib/msynth/plots/genmap.sh -i J_mode_01.txt -z 8 --cbmin -3 --cbmax 3 --cbstep 1.5 --xnodes 360 --ynodes 180 --no-colorbar -t "J@-@~\152@~@- spatial mode 1, 300km" -o J_01_phi_300km.png
#sh /data/palken/mainlib/msynth/plots/genmap.sh -i J_mode_05.txt -z 3 --cbmin -80 --cbmax 80 --cbstep 40 --xnodes 360 --ynodes 180 --no-colorbar -t "J@-r@- spatial mode 5, 300km" -o J_05_r_300km.png
#convert +append J_01_phi_300km.png J_05_r_300km.png out_300.png
