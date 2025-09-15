#!/usr/bin/env bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 outdir"
    exit 1
fi

outdir=$1

genomesdir="data/genomes"
dmddbdir="${outdir}/diamond_db"

mkdir -p $dmddbdir

awk -F ','  'NR > 1 { print $0 }' data/taxid_to_accnum.csv | while IFS=',' read -r taxid accnum; do
    if [[ $accnum == NOTFOUND ]]; then
        continue
    fi
    echo $taxid $accnum
    f=${genomesdir}/${taxid}/ncbi_data/${accnum}/protein.faa
    subdir=${taxid}
    diamond makedb --in "$f" -d "${dmddbdir}/${subdir}" --quiet
done
