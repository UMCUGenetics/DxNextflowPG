#!/usr/bin/env python

from pypgx.sdk import utils as sdk
from sys import argv
import pandas as pd
import argparse

def get_opts():
    p = argparse.ArgumentParser(
        description = "Combine outputs of pypgx/excel calling")
    p.add_argument("--pypgx_dirs",
                   dest='pypgx_dirs',
                   nargs="+")
    p.add_argument("--gene",
                   dest='pgx_gene')
    p.add_argument("--excel_csvs",
                   dest="excel_csvs",
                   nargs="+")

    return p.parse_args()


if __name__ == "__main__":

    args = get_opts()

    """
        Concatenate individual pypgx sample runs with the same PGx gene, in pypgx format
    """
    pypgx_archives = [sdk.Archive.from_file(pypgx_output+"/results.zip").data
                  for pypgx_output in args.pypgx_dirs]
    combined_sample_data = pd.concat(pypgx_archives)
    # Meta data is the same for all samples with the same pgx_gene
    metadata = sdk.Archive.from_file(args.pypgx_dirs[0]+"/results.zip").metadata
    merged_output = sdk.Archive(metadata, combined_sample_data)
    merged_output.to_file(f"{args.pgx_gene}_results.zip")


    """
    Concatenate individual star calling (excel) runs with the same PGx gene
    """
    excel_alleles = pd.concat(
        [pd.read_csv(excel_output, sep="\t", index_col=0)
         for excel_output in args.excel_csvs])

    """
    Concatenate pypgx calls with star calling(excel) calls
    """
    merged = combined_sample_data.join(excel_alleles)
    merged.to_csv(open(f"{args.pgx_gene}.csv", "w"), sep="\t")
