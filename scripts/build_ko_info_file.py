import argparse
import os, sys
import pandas as pd


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--kos_fpath", type=str, 
                        help="file containing list of KOs")
    parser.add_argument("-d", "--datdir", type=str, 
                        help="directory containing KO information files")
    parser.add_argument("-o", "--outfpath", type=str, 
                        help="output filepath")
    return parser.parse_args(args)


def parse_ko_file(filepath):
    entry, symbol, name = None, None, None
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("ENTRY"):
                entry = line.split()[1]  # e.g. "K00360" (ignores trailing "KO")
            elif line.startswith("SYMBOL"):
                symbol = line.split(None, 1)[1].strip()
            elif line.startswith("NAME"):
                name = line.split(None, 1)[1].strip()
    return {"ENTRY": entry, "SYMBOL": symbol, "NAME": name}


def main(args):
    kos_fpath = args.kos_fpath
    datdir = args.datdir
    outfpath = args.outfpath

    
    # Collect list of KOs
    with open(kos_fpath, "r") as f:
        ko_list = [line.strip() for line in f.readlines()]
        ko_filelist = [os.path.join(datdir, ko + ".txt") for ko in ko_list]
        print(ko_list)


    # Collect from all files
    records = [
        parse_ko_file(file) for file in ko_filelist
    ]

    # Put into a dataframe
    df = pd.DataFrame(records)
    df.to_csv(outfpath, sep="\t", index=False)

    print(df.head())
    

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args)
