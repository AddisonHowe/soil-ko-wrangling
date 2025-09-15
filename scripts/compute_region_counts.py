import argparse
import os, sys
import tqdm as tqdm
import numpy as np
import pandas as pd

TAXID_TO_ACCREF_FPATH = "data/taxid_to_accnum.csv"
TAXID_TO_SCAFFOLD_FPATH = "data/taxid_to_scaffold.csv"


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--regions_dir", type=str, required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)
    parser.add_argument("--coverage_dir", type=str, required=True)
    parser.add_argument("--para_coverage_dir", type=str, required=True)
    return parser.parse_args(args)


def main(args):
    outdir = args.outdir
    coverage_dir = args.coverage_dir
    coverage_dir_para = args.para_coverage_dir
    region_tables_dir = args.regions_dir

    # Get list of samples
    sample_filelist = list(
        filter(
            lambda s: s.endswith(".npz"),
            os.listdir(coverage_dir),
        )
    )

    # Get taxids and taxa info
    df_taxids = pd.merge(
        pd.read_csv(TAXID_TO_ACCREF_FPATH, header=0),
        pd.read_csv(TAXID_TO_SCAFFOLD_FPATH, header=0),
    )

    taxid_to_scaffold = {}
    for _, row in df_taxids.iterrows():
        scaffold = row["species"]
        taxid_to_scaffold[row["taxid"]] = scaffold

    # Get list of region tables and construct a map from taxid to the table
    regions_filelist = list(
        filter(lambda s: s.endswith(".tsv"), os.listdir(region_tables_dir))
    )

    taxid_to_region_table = {}
    for fname in regions_filelist:
        taxid = int(fname.removeprefix("identified_regions_").removesuffix(".tsv"))
        fpath = os.path.join(region_tables_dir, fname)
        df = pd.read_csv(fpath, sep="\t", header=0)
        df.insert(0, "taxid", taxid)
        df["seq_length"] = df["stop"] - df["start"]
        taxid_to_region_table[taxid] = df

    # Construct a dataframe for each sample
    print("Constructing dataframes for each sample...")
    os.makedirs(outdir, exist_ok=True)
    for coverage_array in tqdm.tqdm(sample_filelist):
        sample = coverage_array.removeprefix("coverage_arrays_").removesuffix(".npz")
        coverage_arrs = np.load(os.path.join(coverage_dir, coverage_array))
        coverage_arrs_para = np.load(os.path.join(coverage_dir_para, coverage_array))
        df = construct_df_for_sample(
            taxid_to_region_table,
            coverage_arrs,
            coverage_arrs_para,
            taxid_to_scaffold,
        )
        df.to_csv(f"{outdir}/coverage_{sample}.csv", index=False)

    print("Done!")


def is_parabactereoid(taxid):
    return taxid == 2849180


def construct_df_for_sample(
    taxid_to_region_table,
    coverage_arrs,
    coverage_arrs_para,
    taxid_to_scaffold_map,
):
    dfs = []
    for taxid in taxid_to_region_table:
        df0 = taxid_to_region_table[taxid].copy()
        df0["avg_depth"] = -1
        df0["max_depth"] = -1
        if is_parabactereoid(taxid):
            avg_depth = -np.ones(len(df0))
            max_depth = -np.ones(len(df0))
            for i, (idx, row) in enumerate(df0.iterrows()):
                arr_key = row["name"].split(":")[1] + "_Parabacteroides_sp._MSK.9.14"
                base_coverage = coverage_arrs_para[arr_key]
                start, stop = row[["start", "stop"]]
                x = base_coverage[start:stop]
                avg_depth[i] = x.mean()
                max_depth[i] = x.max()
            df0["avg_depth"] = avg_depth
            df0["max_depth"] = max_depth
            dfs.append(df0)

        else:
            arr_key = taxid_to_scaffold_map[taxid]
            avg_depth = -np.ones(len(df0))
            max_depth = -np.ones(len(df0))
            if not pd.isna(arr_key):
                base_coverage = coverage_arrs[arr_key]
                for i, (idx, row) in enumerate(df0.iterrows()):
                    start, stop = row[["start", "stop"]]
                    x = base_coverage[start:stop]
                    avg_depth[i] = x.mean()
                    max_depth[i] = x.max()
                df0["avg_depth"] = avg_depth
                df0["max_depth"] = max_depth
            dfs.append(df0)
    sample_df = pd.concat(dfs)
    sample_df["desc"] = sample_df.pop("desc")
    return sample_df


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args)
