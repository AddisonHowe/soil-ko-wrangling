#!/usr/bin/env bash

genomesdir="data/genomes"
outdir="out/diamond_db"

mkdir -p $outdir

awk -F ','  'NR > 1 { print $0 }' data/taxid_to_accnum.csv | while IFS=',' read -r taxid accnum; do
    if [[ $accnum == NOTFOUND ]]; then
        continue
    fi
    echo $taxid $accnum
    f=${genomesdir}/${taxid}/ncbi_data/${accnum}/protein.faa
    subdir=${taxid}
    diamond makedb --in "$f" -d "${outdir}/${subdir}" --quiet
done
