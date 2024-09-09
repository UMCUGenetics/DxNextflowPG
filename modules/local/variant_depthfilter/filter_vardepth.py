#!/usr/bin/env python3

import argparse
import pandas
from itertools import takewhile

def get_arguments():
    p = argparse.ArgumentParser()
    p.add_argument("--vcf")

    return p.parse_args()


def get_vcf_header(file_path):
    with open(file_path, 'r') as vcf_obj:
        header = list(takewhile(lambda s: s.startwith("#"), vcf_obj))

    return header

if __name__=="__main__":
    args = get_arguments()


    header = get_vcf_header(args.vcf)
    print(header)
