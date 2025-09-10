#!/usr/bin/env bash

indir="out/diamond_res"
outdirbase=out/diamond_hits

mkdir -p $outdirbase
for ko in $(ls ${indir}); do
    echo $ko
    outdir=${outdirbase}/${ko}
    mkdir -p $outdir
    for resfile in ${indir}/$ko/*.tsv; do
        base=$(basename $resfile .tsv)
        taxid=${base/_hits/}
        echo "Finding unique sequence hits for $ko (taxid=$taxid)"
        cut -f2 "$resfile" | sort | uniq > "${outdir}/${taxid}_res_unique.txt"
    done
done
