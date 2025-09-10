#!/usr/bin/env bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 taxid ko"
    exit 1
fi

taxid=$1
ko=$2

queryfile="data/kos/${ko}/${ko}_rep.faa"
dbfile="out/diamond_db/${taxid}.dmnd"
outdir="out/diamond_res/${ko}"
outfname="${taxid}_hits.tsv"

if [[ ! -f $queryfile ]]; then
    echo "Query file ${queryfile} does not exist!"
    exit 1
elif [[ ! -f $dbfile ]]; then
    echo "Diamond database ${dbfile} does not exist!"
    exit 1
fi

mkdir -p $outdir

diamond blastp \
    --query ${queryfile} \
    --db "$dbfile" \
    --out "${outdir}/${outfname}" \
    --outfmt 6 \
    --evalue 1e-5 \
    --max-target-seqs 1 \
    --quiet
