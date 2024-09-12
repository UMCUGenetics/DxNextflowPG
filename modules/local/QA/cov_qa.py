#!/usr/bin/env python3

import pandas as pd
from sys import argv
from os import path


def parse_cov_file(file_path):
     with open(file_path,'r') as file_obj:
         name = path.basename(file_path).replace(".cov","")
         coverage = file_obj.readline().strip()

         return name, coverage


def merge_cov_files(cov_files):
    cov_df = pd.DataFrame(
        [parse_cov_file(path) for path in cov_files],
        columns = ["Sample", "PGx Coverage"])
    return cov_df

if __name__ == "__main__":
    genotypes_df = pd.read_csv(argv[1], sep="\t", index_col = 0)
    cov_df = merge_cov_files(argv[3:])

    merged = pd.merge(
        genotypes_df,
        cov_df,
        how="left",
        left_index=True,
        right_on='Sample')
    merged.set_index("Sample")

    merged.to_csv(argv[2], sep='\t')
