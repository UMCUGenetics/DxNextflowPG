#!/usr/bin/env python3

import argparse
import pandas
import gzip
from itertools import takewhile
import re

info_re = re.compile(r"##INFO=<ID=([A-Z0-9]+)")

def get_arguments():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--vcf")

    p.add_argument(
        "--coverage_field",
        default="DP",
        help="Name of the VCF INFO field to use")

    p.add_argument(
        "--threshold",
        default=15,
        type=int,
        help="Coverage threshold for a SNP to be removed")

    return p.parse_args()


def get_vcf_info_fields(file_path):

    with gzip.open(file_path, 'rt', encoding="utf-8") as vcf_obj:
        headeriter = takewhile(lambda s: s.startswith("#"), vcf_obj)

        for line in headeriter:
            if not line.startswith("##INFO"):
                continue

            match = info_re.search(line)

    yield match.group(1)



if __name__=="__main__":
    args = get_arguments()


    info_fields = parse_vcf_header(args.vcf)

    if "DP"
