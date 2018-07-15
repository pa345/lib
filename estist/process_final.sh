#!/bin/sh
# Create final Est/Ist data file from previously computed files

outfile="$DATAHOME/Indices/ESTIST/Est_Ist_index.pli"

proc_script="$MYLIBHOME/estist/process_year.sh"
datadir="$MYLIBHOME/estist/data"

# Generate new data file for current year
current_year=$(date +%Y)
sh ${proc_script} ${current_year}

echo "Generating ${outfile}"
cat ${datadir}/Est_Ist*.dat | grep -v "\-9999" > ${outfile}
