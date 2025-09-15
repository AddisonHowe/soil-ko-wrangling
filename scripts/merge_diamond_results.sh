#!/usr/bin/env bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 outdir"
    exit 1
fi

outdir=$1

dmddir="${outdir}/diamond_res"
outfpath="${outdir}/dmnd_combined.tsv"

HEADER="taxid\tko\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n"

printf $HEADER > ${outfpath}
for ko in $(ls ${dmddir}); do 
    for f in $(ls ${dmddir}/$ko | grep hits); do 
        taxid=${f/_hits.tsv/}
        cat ${dmddir}/$ko/$f | awk -F'\t' -v ko=$ko -v taxid=$taxid \
            '{print taxid"\t"ko"\t"$0}' >> ${outfpath} 
    done
done
