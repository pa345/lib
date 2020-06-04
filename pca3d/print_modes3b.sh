#!/bin/sh

outdir="modes"

rm -rf ${outdir}
mkdir -p ${outdir}

band="11"

for alt in 110 150 200 250 300 350 400 450; do
  mkdir -p ${outdir}/${alt}
  ./print_modes3b -a ${alt} -J -p ${outdir}/${alt} -f ${band}
done
