#!/usr/bin/env python

from pypgx.sdk import utils as sdk
import pandas as pd



if __name__ == "__main__":

    #output_dirs and gene are substituted by nextflow when running the script as template
    output_dirs = "$pypgx_dirs".split()
    gene = "$pgx_gene"

    # Concatenate individual pypgx sample runs with the same PGx gene, in pypgx format
    pypgx_archives = [sdk.Archive.from_file(pypgx_output+"/results.zip").data
                      for pypgx_output in output_dirs]
    combined_sample_data = pd.concat(pypgx_archives)

    # Output them as csv too
    combined_sample_data.to_csv(open(f"{gene}_mqc.csv", "w"), sep="\t")

    # Meta data is the same for all samples with the same pgx_gene
    metadata = sdk.Archive.from_file(output_dirs[0]+"/results.zip").metadata
    merged_output = sdk.Archive(metadata, combined_sample_data)
    merged_output.to_file(f"{gene}_results.zip")
