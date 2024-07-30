#!/usr/bin/env python


from itertools import chain

# Local
from sample import Sample
import databases



def get_args():
    import argparse

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Clean and format Excel data, and store dbsnp IDstar alleles')

    parser.add_argument(
        '--translation_table',
        dest='excel_path',
        type=str,
        help='Path to the Excel file')

    parser.add_argument(
        "--dbSNP",
        dest='db_snp',
        type=str,
        help='Path to (subsetted) dbSNP vcf')

    parser.add_argument(
        "--vcf",
        dest='vcf_path',
        type=str,
        help="Path to the folder containing the vcf/called zip files")

    parser.add_argument(
        "--pgx_gene",
        dest="pgx_gene",
        type=str,
        help="PGx gene to genotype"
    )

    parser.add_argument(
        "-o",
        dest="outfile",
        type=str,
        help="Path to the output in csv format"
    )
    return parser.parse_args()

def get_lookup_identifier(pair):
    return f"{pair[0]};{pair[1]}"

def get_predominant_genotype(genotypes):
    # todo: make this more sophisticated. It should take into account the phenotype of the genotypes

    wild_type = "*1/*1"
    try:
        genotypes.remove(wild_type)
    except ValueError:
        pass
    
    return genotypes[0]

def has_genotype(genotype, data_subset, rs_gt):

        data_subset_genotype = data_subset[
            data_subset['genotype_show'] == genotype]

        rsID_genotypes = zip(
            data_subset_genotype['variant_id'],
            data_subset_genotype['variant_result'])

        lookup = [get_lookup_identifier(pair) in rs_gt for pair in rsID_genotypes]

        if False in lookup:
            return False
        else:
            return True

def genotype(pgx_gene, data_subset, rs_gt):
    pgx_genotypes = list(set(data_subset['genotype_show']))

    rs_gt = list(chain(*rs_gt))

    
    genotypes = [genotype
                 for genotype in pgx_genotypes if has_genotype(genotype, data_subset, rs_gt)]


    if len(genotypes) == 1:
        return genotypes[0]
    elif len(genotypes) > 1:
        return get_predominant_genotype(genotypes)







def main():
    args = get_args()

    # Clean the data
    excel_conversions = databases.import_clean_excel(args.excel_path)

    

    # Dict to map genomic positions to dbSNP rs ids
    rs_db = databases.import_dbsnp(args.db_snp)

    sample = Sample(args.vcf_path, rs_db)

    with open(args.outfile, 'w') as output_file:

        header = ["Sample", "Gene", "Excel_Genotype"]
        output_file.write("\t".join(header) + '\n')



        pgx_gene = args.pgx_gene
        data_subset = excel_conversions[excel_conversions['gene'] == pgx_gene]


        excel_genotype = genotype(pgx_gene, data_subset, sample.rs_gt)

        outline = [sample.name, pgx_gene, excel_genotype]
        output_file.write("\t".join(outline)+ '\n')
        
    
if __name__ == '__main__':
    main()
