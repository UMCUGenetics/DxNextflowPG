#! venv/bin/python
# standard libraries alphabetic order of main package.
import argparse
import pathlib
import sys
from warnings import warn as warnings_warn

# third party libraries alphabetic order of main package.
import pysam
import vcf as pyvcf
import yaml

# custom libraries alphabetic order
from assets.utils import non_empty_existing_file


def valid_compressed_input_vcf(file):
    non_empty_existing_file(file)
    path_file = pathlib.Path(file)
    if path_file.suffixes != ['.vcf', '.gz'] and path_file.suffix != ".vcf":
        raise OSError("Expected a VCF file (.vcf or .vcf.gz)")
    elif path_file.suffixes == ['.vcf', '.gz']:
        non_empty_existing_file(file)
        non_empty_existing_file(file + ".tbi")
    else:  # .vcf
        pysam.tabix_compress(file, file + ".gz")
        pysam.tabix_index(file + ".gz", preset="vcf")
        file = file + ".gz"
    return file


def valid_translation_table(file):
    if pathlib.Path(file).suffix != ".yaml":
        raise argparse.ArgumentTypeError("Expected a translation table with .yaml extension.")
    return non_empty_existing_file(file)


def parse_arguments_and_check(args_in):
    parser = argparse.ArgumentParser(
        description="Translate variant genotype to a pharmacogentics phenotype.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input", type=valid_compressed_input_vcf, help="File path and name of compressed VCF (.vcf.gz).")
    parser.add_argument("sample", type=str, help="Sample identifier.")
    parser.add_argument("table", type=valid_translation_table, help="File path and name of translation table (.yaml).")
    parser.add_argument(
        "-o", "--output_path", type=non_empty_existing_file, required=False, default=pathlib.Path(__file__).cwd(),
        help="Filepath to store output. (default: %(default)s)"
    )
    parser.add_argument(
        "-p", "--output_prefix", type=str, required=False,
        help="Output prefix to use as output filename. (default: the provided sample identifier)"
    )
    args = parser.parse_args(args_in)

    if not args.output_prefix:
        args.output_prefix = args.sample
    return(args)


def read_vcf(vcf_file):
    try:
        vcf_reader = pyvcf.Reader(filename=vcf_file)
    except StopIteration:
        raise ValueError("File is empty.")
    if not next(vcf_reader, None):
        raise ValueError("File has no records.")
    elif len(next(vcf_reader).samples) != 1:
        raise ValueError("VCF should have (only) one sample.")
    return(vcf_reader)


def read_yaml(translation_file):
    translation_table = yaml.safe_load(translation_file)
    if not translation_table:
        raise ValueError("File is empty.")
    return(translation_table)


def check_required_keys(snp):
    for key in ["chrom", "start", "end", "variant_genotype", "name"]:
        if key not in snp:
            raise ValueError(f"Missing required field {key} for snp in translation table.")


def retrieve_match_all_records(vcf_reader, snp):
    records_match = []
    try:
        records = vcf_reader.fetch(snp.get("chrom"), snp.get("start"), snp.get("end"))
    except ValueError:
        warnings_warn("Remove genotype, fetch has no records for variant {}".format(snp.get("name")))
        return(None)
    for record in records:
        if record.samples[0].gt_bases == snp.get("variant_genotype"):
            records_match.append(True)
        elif record.samples[0]['GT'] == '1/0':
            splitted_genotype = record.samples[0].gt_bases.split("/")[::-1]
            if "/".join(splitted_genotype) == snp.get("variant_genotype"):
                records_match.append(True)
            else:
                records_match.append(False)
        else:
            records_match.append(False)
    if not records_match:  # if list is empty.
        warnings_warn("Remove genotype, fetch has no records for variant {}".format(snp.get("name")))
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
    filter_genotype = []
    genotype_match_per_snp = {}
    for genotype in genotypes:
        if genotype.get("genotype_id") in filter_genotype:
            continue
        for snp in genotype.get("snp"):
            check_required_keys(snp=snp)
            records_match = retrieve_match_all_records(vcf_reader=vcf_reader, snp=snp)
            if not records_match:  # if 'None', remove and filter genotype
                genotype_match_per_snp.pop(genotype.get("genotype_id"), None)
                filter_genotype.append(genotype.get("genotype_id"))
                break
            elif any(records_match):  # when single or all records match
                genotype_match_per_snp.setdefault(genotype.get("genotype_id"), []).append(True)
            elif records_match:  # if none of the records match
                genotype_match_per_snp.setdefault(genotype.get("genotype_id"), []).append(False)
    if not genotype_match_per_snp:  # if all genotypes are removed
        warnings_warn("No genotype match found.")
    return(genotype_match_per_snp)


def write_matched_phenotype(genotype_match_per_snp, translation_table, output_path, output_prefix, sample):
    with open(f"{output_path}/{output_prefix}.txt", 'w') as outfile:
        print("sample\tgenotype_id\tgenotype_match\tphenotype_id\tphenotype_match", file=outfile)
        for genotype_id, genotype_match in genotype_match_per_snp.items():
            if all(genotype_match):
                print("{sample}\t{gt_id}\t{gt}\t{pt_nr}\t{pt}".format(
                    sample=sample,
                    gt_id=genotype_id,
                    gt=translation_table[genotype_id].get("gene_genotype"),
                    pt_nr=translation_table[genotype_id].get("phenotype_id"),
                    pt=translation_table[genotype_id].get("phenotype_name"),
                ), file=outfile)


def main(translation_file, vcf_file, output_path, output_prefix, sample):
    vcf_reader = read_vcf(vcf_file)
    translation_table = read_yaml(translation_file)
    translation_file.close()
    genotype_match_per_snp = retrieve_match_snp_genotype(vcf_reader=vcf_reader, genotypes=translation_table.values())
    if genotype_match_per_snp:
        write_matched_phenotype(
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
