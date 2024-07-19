#!/usr/bin/env python3

import gzip
import random
import os

from sys import argv

def load_vcf_file(vcf_file):
    if vcf_file.endswith('.gz'):
        return gzip.open(vcf_file, 'rt')
    else:
        return open(vcf_file, 'r')

def get_samples(vcf_file):
    with load_vcf_file(vcf_file) as file:
        for line in file:
            if line.startswith('#CHROM'):
                return line.strip().split('\t')[9:]

def randomly_select_samples(samples, num_samples=10):
    return random.sample(samples, num_samples)

def write_output_vcf(vcf_file, selected_sample, output_file):
    with load_vcf_file(vcf_file) as file, open(output_file, 'w') as out_vcf:
        for line in file:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    headers = line.strip().split('\t')
                    out_vcf.write('\t'.join(headers[:9] + [selected_sample]) + '\n')
                else:
                    out_vcf.write(line)
            else:
                columns = line.strip().split('\t')
                sample_index = headers.index(selected_sample)
                genotype = columns[sample_index].split(':')[0]  # Extract the GT field

                # Mock values for AD, DP, and AF fields
                ad_value = '10,5'
                dp_value = '15'
                af_value = '0.33'
                genotype_str = f"{genotype}:{ad_value}:{dp_value}:{af_value}"

                out_vcf.write('\t'.join(columns[:7] + ['Phased'] + ['GT:AD:DP:AF'] + [genotype_str]) + '\n')

def main(vcf_file, output_dir, num_samples=10):
    samples = get_samples(vcf_file)
    selected_samples = randomly_select_samples(samples, num_samples)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for sample in selected_samples:
        output_file = os.path.join(output_dir, f"{sample}_genotypes.vcf")
        write_output_vcf(vcf_file, sample, output_file)

# Example usage:
vcf_file = argv[1]  # 'variants_regions.vcf.gz'
output_dir = argv[2] # 'output_genotypes/'

main(vcf_file, output_dir)
