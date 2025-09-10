#!/usr/bin/env bash

outdir="out/identified_regions_best_match"
diamond_res_fpath="out/dmnd_combined_best_match.tsv"

awk -F ','  'NR > 1 { print $0 }' data/taxid_to_accnum.csv | while IFS=',' read -r taxid accnum; do
    echo $taxid $accnum
    prefix=${dst}/${taxid}/ncbi_data
    if [[ $accnum == NOTFOUND ]]; then
        continue
    fi
    python scripts/locate_regions_in_genome.py \
        -t ${taxid} \
        -a data/genomes/${taxid}/ncbi_data/${accnum}/genomic.gff \
        -d ${diamond_res_fpath} \
        -o ${outdir}
done

echo Done!
