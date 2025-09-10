#!/bin/bash

datdir="data/genomes"

cdfname="cds_from_genomic.fna"
outfname="cd_names.txt"

awk -F ','  'NR > 1 { print $0 }' taxid_to_accnum.csv | while IFS=',' read -r taxid accnum; do
    echo $taxid $accnum
    prefix=${datdir}/${taxid}/ncbi_data
    if [[ $accnum == NOTFOUND ]]; then
        continue
    else
        cat ${prefix}/${accnum}/${cdfname} | grep ">" > ${datdir}/${taxid}/${outfname}
    fi
done


cdfname="protein.faa"
outfname="protein_names.txt"

awk -F ','  'NR > 1 { print $0 }' taxid_to_accnum.csv | while IFS=',' read -r taxid accnum; do
    echo $taxid $accnum
    prefix=${datdir}/${taxid}/ncbi_data
    if [[ $accnum == NOTFOUND ]]; then
        continue
    else
        cat ${prefix}/${accnum}/${cdfname} | grep ">" > ${datdir}/${taxid}/${outfname}
    fi
done
