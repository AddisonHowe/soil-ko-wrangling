#!/usr/bin/env bash

outdirbase=out/diamond_res_uniq
mkdir -p $outdirbase
for ko in $(ls out/diamond_res); do
    echo $ko
    outdir=${outdirbase}/${ko}
    mkdir -p $outdir
    for resfile in out/diamond_res/$ko/*.tsv; do
        base=$(basename $resfile .tsv)
        taxid=${base/_hits/}
        echo "Finding unique sequence hits for $ko (taxid=$taxid)"
        cut -f2 "$resfile" | sort | uniq > "${outdir}/${taxid}_res_unique.txt"
    done
done
