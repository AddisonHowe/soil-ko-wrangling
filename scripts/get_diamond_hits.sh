#!/usr/bin/env bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 outdir"
    exit 1
fi

outdir=$1

dmdresdir="${outdir}/diamond_res"
dmdhitdir="${outdir}/diamond_hits"

mkdir -p $dmdhitdir
for ko in $(ls ${dmdresdir}); do
    echo $ko
    subdir=${dmdhitdir}/${ko}
    mkdir -p $subdir
    for resfile in ${dmdresdir}/$ko/*.tsv; do
        base=$(basename $resfile .tsv)
        taxid=${base/_hits/}
        echo "Finding unique sequence hits for $ko (taxid=$taxid)"
        cut -f2 "$resfile" | sort | uniq > "${subdir}/${taxid}_res_unique.txt"
    done
done
