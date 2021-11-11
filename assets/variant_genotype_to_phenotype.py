#! venv/bin/python
# standard libraries alphabetic order of main package.
import argparse
from errno import ENOENT as errno_ENOENT
from os import strerror as os_strerror
import pathlib
import sys
from warnings import warn as warnings_warn

# third party libraries alphabetic order of main package.
import pysam
import vcf as pyvcf
import yaml

def parse_arguments_and_check(args_in):
    parser = argparse.ArgumentParser(description="Translate variant genotype to a pharmacogentics phenotype.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", type=str, help="File path and name of compressed VCF (.vcf.gz).")
    parser.add_argument("sample", type=str, help="Sample identifier.")
    parser.add_argument("table", type=str, help="File path and name of translation table (.yaml).")
    parser.add_argument("-o", "--output_path", type=str, required=False, default=os.getcwd(), help="File path to store output.")
    parser.add_argument("-p", "--output_prefix", type=str, required=False,
                        help="Output prefix to use as output filename. (default: the provided sample identifier)")
    args = parser.parse_args(args_in)
    if not args.output_prefix:
        args.output_prefix = args.sample
    if not args.table.endswith(".yaml"):
        raise TypeError()
    if not args.input.endswith(".vcf.gz"):
        pysam.tabix_compress(args.input, args.input + ".gz")
        pysam.tabix_index(args.input + ".gz", preset="vcf")
        args.input = args.input + ".gz"
    for input_file_or_dir in [args.table, args.input, args.input+".tbi", args.output_path]:
        if not os.path.isfile(input_file_or_dir) and not os.path.isdir(input_file_or_dir):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), input_file_or_dir)
    return(args)


def read_vcf(vcf_file):
    try:
        vcf_reader = pyvcf.Reader(filename=vcf_file)
    except StopIteration:
        raise ValueError("File is empty.")
    if not next(vcf_reader, None):
        raise ValueError("File has no records.")
    return(vcf_reader)


def read_table(translation_file):
    with open(translation_file) as translation_file:
        translation_table = yaml.safe_load(translation_file)
    if not translation_table:
        raise ValueError("File is empty.")
    return(translation_table)
    

def check_required_keys(snp):
    for key in ["chrom", "start", "end", "variant_genotype", "name"]:
        if key not in snp:
            raise ValueError("Missing required field {} for snp in translation table.".format(key))


def retrieve_match_all_records(vcf_reader, snp):
    records_match = []
    try:
        records = vcf_reader.fetch(snp.get("chrom"), snp.get("start"), snp.get("end"))
    except ValueError:
        warnings_warn("Remove genotype, fetch has no records for variant {}".format(snp.get("name")))
        return(None)
    for record in records:
        if record.samples[0].gt_bases == snp.get("variant_genotype") :
            records_match.append(True)
        elif record.samples[0]['GT'] == '1/0':
            splitted_genotype = record.samples[0].gt_bases.split("/")[::-1]
            if "/".join(splitted_genotype) == snp.get("variant_genotype"):
                records_match.append(True)
            else:
                records_match.append(False)
        else:
            records_match.append(False)
    if not records_match:
        print("Remove genotype, fetch has no records for variant {}".format(snp.get("name")))
        return(None)
    return(records_match)


def retrieve_match_snp_genotype(vcf_reader, genotypes):
    '''
    A genotype/phenotype is represented by at least one site of interest, and often multiple.
    A single site of interest can map to multiple probes on the array.
    For each genomic location of a site of interest, are VCF records retrieved.
    At least one VCF record should match the required snp genotype to conclude the sample
    matches.
    '''
    genotype_match_per_snp = {}
    for genotype in genotypes:
        for snp in genotype.get("snp"):
            check_required_keys(snp=snp)
            records_match = retrieve_match_all_records(vcf_reader=vcf_reader, snp=snp)
            if not records_match:
                genotype_match_per_snp.pop(genotype.get("genotype_id"), None)
                break # TODO: refactor code to break nested loop.
            elif any(records_match):
                genotype_match_per_snp.setdefault(genotype.get("genotype_id"), []).append(True)
            elif records_match:
                genotype_match_per_snp.setdefault(genotype.get("genotype_id"), []).append(False)
                
    return(genotype_match_per_snp)


def retrieve_match_phenotype_and_write(genotype_match_per_snp, translation_table, output_path, output_prefix, sample):
    with open("{path}/{prefix}.txt".format(path=output_path, prefix=output_prefix), 'w') as outfile:
        print("sample\tgenotype_match\tphenotype_tsitnr\tphenotype_match", file=outfile)
        for genotype_id, genotype_match in genotype_match_per_snp.items():
            if all(genotype_match):
                print("{sample}\t{gt}\t{pt_nr}\t{pt}".format(
                    sample=sample,
                    gt=translation_table[genotype_id].get("genotype_autogenerated_name"),
                    pt_nr=translation_table[genotype_id].get("gstandard_tsitnr"),
                    pt=translation_table[genotype_id].get("gstandard_phenotype_name_thnm50"),
                ), file=outfile)


def main(translation_file, vcf_file, output_path, output_prefix, sample):
    vcf_reader = read_vcf(vcf_file)
    translation_table = read_table(translation_file)
    genotype_match_per_snp = retrieve_match_snp_genotype(vcf_reader=vcf_reader, genotypes=translation_table.values())
    retrieve_match_phenotype_and_write(
        genotype_match_per_snp=genotype_match_per_snp,
        translation_table=translation_table,
        output_path=output_path,
        output_prefix=output_prefix,
        sample=sample,
        )


if __name__ == '__main__':
    args = parse_arguments_and_check(args_in=sys.argv[1:])
    main(
        translation_file=args.table,
        vcf_file=args.input,
        output_path=args.output_path,
        output_prefix=args.output_prefix,
        sample=args.sample,
        )
