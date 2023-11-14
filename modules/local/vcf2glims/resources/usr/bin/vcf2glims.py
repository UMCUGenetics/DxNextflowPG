#! /usr/bin/env python

import argparse

import vcfpy


def vcf2csv(args):
    reader = vcfpy.Reader.from_path(args.vcf_file)

    # print header
    print('sample', 'sequencing_run', 'chrom', 'pos', 'id', 'genotype', sep=',')

    for record in reader:
        sample_call = record.calls[0]  # Assume single sample vcf

        # Set genotype separator
        gt_sep = '/'
        if sample_call.is_phased:
            gt_sep = '|'

        # create genotype using genotype separator
        sample_gt = gt_sep.join(sample_call.gt_bases)

        if args.min_dp and sample_call.data['DP'] < args.min_dp:
            sample_gt = 'qc_fail_dp'
        elif args.min_gq and sample_call.is_variant and sample_call.data['GQ'] < args.min_gq:
            sample_gt = 'qc_fail_gq'
        elif args.min_rgq and not sample_call.is_variant and sample_call.data['RGQ'] < args.min_rgq:
            sample_gt = 'qc_fail_rgq'

        print(
            sample_call.sample,
            args.sequencing_run,
            record.CHROM,
            record.POS,
            record.ID[0] if record.ID else '',  # print record.ID unless empty
            sample_gt,
            sep=','
        )


if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser(description="")
    parser.set_defaults(func=vcf2csv)
    # Required arguments
    parser.add_argument("sequencing_run", help="")
    parser.add_argument("vcf_file", type=str, help="")
    # Optional arguments
    parser.add_argument("--min_dp", type=int, help="Minimum DP")
    parser.add_argument("--min_gq", type=int, help="Minimum GQ")
    parser.add_argument("--min_rgq", type=int, help="Minimum RGQ")

    args = parser.parse_args()
    args.func(args)