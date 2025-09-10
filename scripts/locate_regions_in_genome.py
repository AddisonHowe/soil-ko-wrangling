import argparse
import os, sys
import csv
import numpy as np
import pandas as pd


ANNOTATION_DESC_COL = 8
ANNOTATION_KEY_OF_INTEREST = "Name"


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--taxid", type=int, required=True)
    parser.add_argument("-a", "--annotation_fpath", type=str, required=True)
    parser.add_argument("-d", "--diamond_fpath", type=str, required=True,
                        help="combined diamond results filtered to unique gene ids")
    parser.add_argument("-o", "--outdir", type=str, required=True)
    parser.add_argument("-v", "--verbosity", type=int, default=1)
    return parser.parse_args(args)


def main(args):

    taxid = args.taxid
    annotation_fpath = args.annotation_fpath
    diamond_fpath = args.diamond_fpath
    outdir = args.outdir
    verbosity = args.verbosity

    if verbosity:
        print(f"TAXID: {taxid}")
        print(f"diamond results: {diamond_fpath}")

    df = pd.read_csv(
        diamond_fpath,
        sep="\t"
    )


    # Subset data to taxid of interest.
    df = df[df["taxid"] == taxid]

    if verbosity > 1:
        print(df)

    # Build map from KO to regions ids
    ko_to_region_ids = {}
    for ko in df["ko"].unique():
        if verbosity > 1:
            print(ko)
        ko_to_region_ids[ko] = df[df["ko"] == ko]["sseqid"].unique()
    if verbosity > 1:
        print(ko_to_region_ids)
    os.makedirs(outdir, exist_ok=True)
    
    identified_regions_by_ko = {ko: [] for ko in ko_to_region_ids}
    with open(annotation_fpath, "r") as annotations:
        csvreader = csv.reader(annotations, delimiter="\t")    
        for row in csvreader:
            # Skip comments and non CDS rows.
            if row[0][0] == "#" or row[2] != "CDS":
                continue
            # Look at description column, which contains a number of pairs of 
            # the form key=value, separated by semicolons.
            desc_str = row[ANNOTATION_DESC_COL]
            desc_items = desc_str.split(";")
            d = {k: v for k, v in [item.split("=") for item in desc_items]}
            # Get the key `name` and check if it's one of the sequence ids
            if ANNOTATION_KEY_OF_INTEREST not in d:
                continue
            name = d[ANNOTATION_KEY_OF_INTEREST]
            for ko in ko_to_region_ids:
                if name in ko_to_region_ids[ko]:
                    start = row[3]
                    stop = row[4]
                    reg = (ko, name, start, stop, desc_str)
                    identified_regions_by_ko[ko].append(reg)
                    
    with open(f"{outdir}/identified_regions_{taxid}.tsv", "w") as fout:
        header = "\t".join(["ko", "name", "start", "stop", "desc"])
        fout.write(header + "\n")
        for ko in identified_regions_by_ko:
            for reg in identified_regions_by_ko[ko]:
                fout.write("\t".join(reg) + "\n")    

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args)
