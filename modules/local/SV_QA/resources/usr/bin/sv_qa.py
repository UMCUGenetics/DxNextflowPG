#!/usr/bin/env python3

"""Performs quality control on SV calls by comparing the specific call to previously calculated frequencies. Rare variants will receive a warning"""

from sys import argv
from glob import glob
from os import path
import pandas as pd


def get_PGx_genes_frequency_dict(frequency_dir):
    """Create a named dictionary containing pharmacogenes with calculated frequencies and their file paths. e.g., {CYP2D6: cyp2d6_frequencies.csv}

    """

    return {path.basename(x).split("_")[0]: x for x in glob(frequency_dir + "/*_freqs.csv")}


def get_frequency_dict(freq_files, PGx_gene, sep=";"):
    """Turn the frequency file(s) into a dictionary with format {SV_name:frequency}
    """
    freq_path = freq_files[PGx_gene]
    freq_df = pd.read_csv(freq_path, sep=sep)

    return freq_df.set_index(PGx_gene).T.to_dict("list")


def get_annotation(row, frequencies, freq_threshold=5):
    """Lookup the SV call for a row in the frequencies dictionary, and give a warning for rare SVs"""
    cnv = row["CNV"]

    try:
        count, freq = frequencies[cnv]
        freq = float(freq.replace(",", "."))

    # If a SV is not in the frequency call I would be very surprised, but if it happens it willt return a freq of 0.0 (and therefore give a warning)
    except KeyError:
        freq = 0.0

    if freq < freq_threshold:
        annotation = "Rare CNV"
    else:
        annotation = "N/A"
    row["CNV Annotation Warnin"] = annotation
    return row


if __name__ == "__main__":
    # Parse the frequency directory for frequency files.
    frequency_files = get_PGx_genes_frequency_dict(argv[1])

    # Output csv from pypgx, which includes the SV calls in a separate column
    genotypes_df = pd.read_csv(argv[2], sep="\t")

    # Extract pharmacogene name from the pypgx output filename
    PGx_gene = path.basename(argv[2]).replace(".csv", "").replace("_mqc", "")

    try:
        # Parse the PGx_gene frequency file for precalculated SV fruquencies
        freq_dict = get_frequency_dict(frequency_files, PGx_gene)

        # Row-wise apply on the pypgx output (each row can contain a different sample and thus a different SV call)
        # get_annotation adds a new column with a SV annotation called 'annotation'
        genotypes_df = genotypes_df.apply(get_annotation, axis=1, frequencies=freq_dict)
    except KeyError:
        # When the frequency table is not available for the PGx gene
        pass

    # When the frequency table is unavailable the original df is outputted
    genotypes_df.to_csv(argv[3], sep="\t")
