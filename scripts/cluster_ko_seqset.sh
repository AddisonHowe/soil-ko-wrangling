#!/usr/bin/env bash


if [ "$#" -ne 1 ]; then
    echo "Usage: $0 ko"
    exit 1
fi

ko=$1

outdir=data/kos/$ko

# Cluster sequences to get representatives
infile="${outdir}/${ko}.faa"
repfile="${outdir}/${ko}_rep.faa"
echo "Clustering sequences with CD-HIT..."
cd-hit -i "$infile" -o "$repfile" -c 0.9 -n 5 -d 0
