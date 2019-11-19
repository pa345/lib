#!/bin/sh
#
# Run stage1 for Cryosat and all years

for year in `seq 2014 2018`; do
  screen -d -m -S p${year} sh stage1_cryosat.sh ${year}
done
