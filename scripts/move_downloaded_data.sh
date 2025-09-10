#!/bin/bash

src="downloads"   # where your TAXID folders are
dst="data/genomes"     # target root directory
mkdir -p "$dst"

awk -F ','  'NR > 1 { print $0 }' taxid_to_accnum.csv | while IFS=',' read -r taxid accnum; do
    echo $taxid $accnum
    prefix=${dst}/${taxid}/ncbi_data
    mkdir -p $prefix
    if [[ $accnum == NOTFOUND ]]; then
        continue
    else
        mv ${src}/${taxid}/data/${accnum} ${prefix}/${accnum}
    fi
done
