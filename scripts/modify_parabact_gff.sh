#!/usr/bin/env bash


genome_annot_fpath=data/genomes/2849180/ncbi_data/GCF_019041495.1/genomic.gff
protein_to_scaffold_fpath=data/genomes/2849180/protein_to_scaffold.tsv
protein_fpath_orig=data/genomes/2849180/ncbi_data/GCF_019041495.1/protein.faa
protein_fpath_new=data/genomes/2849180/ncbi_data/GCF_019041495.1/protein_new.faa

awk -F'\t' '$3=="CDS" {
    pid=$9; scaf=$1;
    sub(/.*protein_id=/,"",pid); sub(/;.*/,"",pid);
    sub(/ .*/,"",scaf);
    print pid "\t" scaf;
}' ${genome_annot_fpath} > ${protein_to_scaffold_fpath}


awk 'NR==FNR {map[$1]=$2; next}
    /^>/ {
        # remove ">" and split header
        hdr=$0; sub(/^>/,"",hdr); split(hdr,a," ");
        pid=a[1];
        if (pid in map) {
            print ">" pid ":" map[pid];
        } else {
            print ">" pid;
        }
        next
    }
    {print}' ${protein_to_scaffold_fpath} ${protein_fpath_orig} > ${protein_fpath_new}
