import os, sys
import time
import numpy as np
import argparse
import tqdm as tqdm
from Bio.KEGG import REST


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", "--ko", type=str, required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)
    parser.add_argument("-v", "--verbosity", type=int, default=1)
    parser.add_argument("-b", "--batch_size", type=int, default=10)
    parser.add_argument("--sleep", type=float, default=0.01)
    parser.add_argument("--pbar", action="store_true")
    return parser.parse_args(args)


def main(args):
    # Step 1: Get all genes linked to a KO
    ko = args.ko
    outdir = args.outdir
    batch_size = args.batch_size
    sleep_time = args.sleep
    verbosity = args.verbosity
    use_pbar = args.pbar

    link_result = REST.kegg_link("genes", ko).read()

    # Parse gene IDs (right-hand side of tab-delimited lines)
    gene_ids = [line.split("\t")[1] for line in link_result.strip().split("\n")]
    print(f"Found {len(gene_ids)} genes for {ko}")

    # Step 2: Fetch amino acid sequences in FASTA format (in batches)
    faa_sequences = fetch_faa(
        gene_ids, 
        batch_size=batch_size, 
        sleep_time=sleep_time,
        verbosity=verbosity,
        pbar=use_pbar
    )

    # Step 3: Save to file
    os.makedirs(outdir, exist_ok=True)
    outfpath = f"{outdir}/{ko}.faa"
    with open(outfpath, "w") as f:
        f.write(faa_sequences)

    print(f"Saved {len(gene_ids)} sequences to {outfpath}")


def fetch_faa(gene_list, batch_size=10, sleep_time=0, verbosity=1, pbar=True):
        if batch_size > 10:
             raise RuntimeError("Max batch size for KEGG API is 10")
        sequences = []
        for i in tqdm.trange(0, len(gene_list), batch_size, disable=not pbar):
            if (not pbar) and verbosity:
                 print(f"Fetching batch {i+1}/{len(gene_list)}")
            batch = gene_list[i:i+batch_size]
            batch_str = "+".join(batch)
            time.sleep(sleep_time)
            faa = REST.kegg_get(batch_str, "aaseq").read()
            sequences.append(faa)
        return "".join(sequences)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args)
