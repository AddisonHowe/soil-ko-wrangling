#!/usr/bin/env bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 outdir"
    exit 1
fi

outdir=$1

outdirbase="${outdir}/identified_regions"
diamond_res_fpath="${outdir}/dmnd_combined_top_hits.tsv"

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
        -o ${outdirbase}
done

echo Done!
