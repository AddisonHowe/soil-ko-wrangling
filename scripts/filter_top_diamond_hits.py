import argparse
import sys
import pandas as pd

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infpath", type=str, required=True)
    parser.add_argument("-o", "--outfpath", type=str, required=True)
    return parser.parse_args(args)


def find_top_hit(df):
        assert len(df["sseqid"].unique()) == 1, "df should have only one sseqid"
        # Determine the minimal evalue
        min_e = df["evalue"].min()
        # Filter rows with that minimum evalue
        tempdf = df[df["evalue"] == min_e]
        # Select the one with the maximum pident
        best_idx = tempdf["pident"].idxmax()
        best_ko = tempdf.loc[best_idx, "ko"]
        best_pident = tempdf.loc[best_idx, "pident"]
        return best_idx, best_ko, min_e, best_pident


def main(args):

    dmd_res_combined_fpath = args.infpath
    outfpath = args.outfpath

    DF = pd.read_csv(
        dmd_res_combined_fpath,
        sep="\t"
    )

    df_list = []
    for taxid in DF["taxid"].unique():
        taxdf = DF[DF["taxid"] == taxid]
        ko_counts_by_geneid = taxdf.groupby("sseqid")["ko"].nunique()
        uniques = ko_counts_by_geneid[ko_counts_by_geneid == 1]
        nonuniques = ko_counts_by_geneid[ko_counts_by_geneid > 1]
        # Process hits associated to a single KO
        for sseq in uniques.index:
            subdf = taxdf[taxdf["sseqid"] == sseq]
            best_idx, best_ko, best_evalue, best_pident = find_top_hit(subdf)
            df_list.append(subdf.loc[best_idx])

        # Process hits associated to multiple KOs
        for sseq in nonuniques.index:
            subdf = taxdf[taxdf["sseqid"] == sseq]
            evalue_min = subdf.groupby("ko")["evalue"].min()
            df_minevalue = subdf[subdf["evalue"].isin(evalue_min)]
            pident_max = df_minevalue.groupby("ko")["pident"].max()
            df_filtered = df_minevalue[df_minevalue["pident"].isin(pident_max)]
            best_idx, _, _, _ = find_top_hit(df_filtered)
            df_list.append(df_filtered.loc[best_idx])

    dfout = pd.DataFrame(df_list)

    dfout.to_csv(
        outfpath,
        index=None,
        sep="\t",
    )

    

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args)
