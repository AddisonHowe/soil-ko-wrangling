#!/usr/bin/env bash

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 taxid ko outdir evalue_thresh pident_thresh"
    exit 1
fi

taxid=$1
ko=$2
outdir=$3
evalue_thresh=$4  # 1e-5
pident_thresh=$5  # 50

queryfile="data/kos/${ko}/${ko}_rep.faa"
dbfile="${outdir}/diamond_db/${taxid}.dmnd"
outdirbase="${outdir}/diamond_res/${ko}"
outfname="${taxid}_hits.tsv"

if [[ ! -f $queryfile ]]; then
    echo "Query file ${queryfile} does not exist!"
    exit 1
elif [[ ! -f $dbfile ]]; then
    echo "Diamond database ${dbfile} does not exist!"
    exit 1
fi

mkdir -p $outdirbase

diamond blastp \
    --query ${queryfile} \
    --db "$dbfile" \
    --out "${outdirbase}/${outfname}" \
    --outfmt 6 \
    --evalue ${evalue_thresh} \
    --id ${pident_thresh} \
    --max-target-seqs 1 \
    --quiet
