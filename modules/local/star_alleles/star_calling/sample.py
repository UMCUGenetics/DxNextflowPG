#!/usr/bin/env python3

import os
import gzip

from databases import get_vcf_identifier

class Sample:

    def __init__(self, vcf, rs_db):
        self.vcf = vcf
        self.rs_gt = []
        self.pypgx_star_alleles = None
        self.directory = os.path.dirname(vcf)
        self.basename = os.path.basename(vcf).replace(".vcf.gz", "")
        self.name = self.basename

        for data_line in self._parse_vcf():
            rs_ID_genotypes = self.import_rsIDs(data_line, rs_db)
            self.rs_gt += rs_ID_genotypes

        self.rs_gt = set(self.rs_gt)

    def _parse_vcf(self):
        with gzip.open(self.vcf) as vcf_file:
            for line in vcf_file:
                line = line.decode()
                if line.startswith("#"):
                    continue
                else:
                    yield line.strip().split()

    @staticmethod
    def import_rsIDs(data_line, rs_db):
        identifier = get_vcf_identifier(data_line)

        try:
            rs_ids = rs_db[identifier]
        except KeyError:
            return []

        genotype = data_line[8].split(":")[0]

        for rs_id in rs_ids:
            rs_ID_genotypes = []
            try:
                genotype_nuc = GT_to_Nucleotides((data_line[3], data_line[4]), genotype)
                rsID_genotype = f"{rs_id};{genotype_nuc}"
                # logging.debug(rsID_genotype)

                rs_ID_genotypes.append(rsID_genotype)
            except IndexError: # chrX has only one haplotype
                pass
        return rs_ID_genotypes

    # def _import_pypgx_star_alleles(self):
    #     zip_dir = f"{self.directory}/{self.basename}_*star_alleles.zip"
    #     zipfile = glob(zip_dir)[0]
    #     df = extract_tsv_from_zip(zipfile)
    #     self.pypgx_star_alleles = df

    # def get_pypgx_star_allele(self):
    #     df = self.pypgx_star_alleles
    #     hap1 = df['Haplotype1'].iloc[0]
    #     hap2 = df["Haplotype2"].iloc[0]
    #     return f"{hap1}/{hap2}"


# def extract_tsv_from_zip(zip_file):
#     with zipfile.ZipFile(zip_file, 'r') as zip_ref:
#         # Zoek naar de data.tsv in een submap
#         for file_info in zip_ref.infolist():
#             if file_info.filename.endswith('data.tsv'):
#                 with zip_ref.open(file_info) as file:
#                     # Lees het bestand in een pandas DataFrame
#                     df = pd.read_csv(file, delimiter='\t')
#                     return df
#     return None

def GT_to_Nucleotides(alleles, GT):
    """ Convert the phased genotype (e.g., 0|0) to the corresponding nucleotides
    in the format of the PGx tanslation table/excel file (e.g. A:A)
    """
    if "|" in GT:
        GT_sep = "|"
    else:
        GT_sep = "/"
    gt_s = [int(x) for x in GT.split(GT_sep)]

    allele_a = alleles[gt_s[0]]
    allele_b = alleles[gt_s[1]]
    return f"{allele_a}:{allele_b}"
