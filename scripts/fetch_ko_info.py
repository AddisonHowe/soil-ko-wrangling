import os, sys
import time
import argparse
import tqdm as tqdm
from Bio.KEGG import REST


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", type=str, required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)
    parser.add_argument("-v", "--verbosity", type=int, default=1)
    parser.add_argument("--sleep", type=float, default=0.5)
    parser.add_argument("--pbar", action="store_true")
    return parser.parse_args(args)


def main(args):
    infile = args.infile
    outdir = args.outdir
    verbosity = args.verbosity
    sleep = args.sleep
    use_pbar = args.pbar
    
    with open(infile, "r") as f:
        lines = f.readlines()
        ko_list = [line.strip() for line in lines]
    
    os.makedirs(outdir, exist_ok=True)

    results = {}
    for i in tqdm.trange(len(ko_list), disable=not use_pbar):
        ko = ko_list[i]
        if (not use_pbar) and verbosity:
            print(ko)
        results[ko] = REST.kegg_get(f"ko:{ko}").read()
        time.sleep(sleep)
        outfpath = os.path.join(outdir, ko + ".txt")
        with open(outfpath, "w") as f:
            f.write(results[ko])


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args)
