#!/usr/bin/env python3

from sys import argv
from glob import glob
from os import path
import pandas as pd


def get_PGx_genes_frequency_dict(frequency_dir):
    return {
        path.basename(x).split("_")[0]: x for x in glob(frequency_dir + "/*_freqs.csv")
    }


def get_frequency_dict(freq_files, PGx_gene, sep=";"):
    freq_path = freq_files[PGx_gene]
    freq_df = pd.read_csv(freq_path, sep=sep)

    return freq_df.set_index(PGx_gene).T.to_dict("list")


def get_annotation(row, frequencies, freq_threshold=1):
    cnv = row["CNV"]

    try:
        count, freq = frequencies[cnv]
        freq = float(freq.replace(",", "."))
    except KeyError:
        freq = 0.0

    if freq < freq_threshold:
        annotation = "Warning: rare CNV detected"
    else:
        annotation = "N/A"
    row["Warning"] = annotation
    return row


if __name__ == "__main__":
    frequency_files = get_PGx_genes_frequency_dict(argv[1])

    genotypes_df = pd.read_csv(argv[2], sep="\t")
    PGx_gene = path.basename(argv[2]).replace(".csv", "")

    try:
        freq_dict = get_frequency_dict(frequency_files, PGx_gene)
        genotypes_df = genotypes_df.apply(get_annotation, axis=1, frequencies=freq_dict)
    except KeyError:
        # When the frequency table is not available for the PGx gene
        pass

    # When the frequency table is unavailable the original df is outputted
    genotypes_df.to_csv(argv[3], sep="\t")
