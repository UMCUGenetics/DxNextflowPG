import os
import pytest
from pytest_reqs import check_requirements
import vcf as pyvcf

from assets.variant_genotype_to_phenotype import retrieve_match_all_records
@pytest.fixture(scope="module", autouse=True)
def get_vcf_reader(setup_and_get_test_path):
    vcf_reader = pyvcf.Reader(filename=setup_and_get_test_path + "/fake_000000000000_R00C00.vcf.gz")
    return(vcf_reader)


VCF_PAIR = (os.path.realpath("./tests/test_data/fake_000000000000_R00C00.vcf.gz"),
              os.path.realpath("./tests/test_data/fake_000000000000_R00C00.vcf.gz.tbi"))


@pytest.fixture(scope="module")
def get_vcf_reader():
    vcf = VCF_PAIR[0]
    vcf_reader = pyvcf.Reader(filename=vcf)
    return(vcf_reader)


class TestSnpGenotypes():
    def retrieve_match_and_assert(self, vcf_reader, snp_dir, exp_length, exp_bool):
        matches = retrieve_match_all_records(vcf_reader=vcf_reader, snp=snp_dir)
        assert(len(matches) == exp_length)
        assert(all(matches) == exp_bool)

    def test_single_snp_ref_ref_match(self, get_vcf_reader):
        snp_dir = {"chrom": "1", "start": 10000000, "end": 10000001, "variant_genotype": "C/C"}
        self.retrieve_match_and_assert(vcf_reader=get_vcf_reader, snp_dir=snp_dir, exp_length=1, exp_bool=True)

    def test_single_snp_ref_alt_match(self, get_vcf_reader):
        snp_dir = {"chrom": "1", "start": 10000001, "end": 10000002, "variant_genotype": "C/T"}
        self.retrieve_match_and_assert(vcf_reader=get_vcf_reader, snp_dir=snp_dir, exp_length=1, exp_bool=True)

    def test_single_snp_alt_ref_match(self, get_vcf_reader):
        snp_dir = {"chrom": "1", "start": 10000002, "end": 10000003, "variant_genotype": "C/T"}
        self.retrieve_match_and_assert(vcf_reader=get_vcf_reader, snp_dir=snp_dir, exp_length=1, exp_bool=True)

    def test_single_snp_alt_alt_match(self, get_vcf_reader):
        snp_dir = {"chrom": "1", "start": 10000003, "end": 10000004, "variant_genotype": "T/T"}
        self.retrieve_match_and_assert(vcf_reader=get_vcf_reader, snp_dir=snp_dir, exp_length=1, exp_bool=True)

    def test_single_snp_mismatch(self, get_vcf_reader):
        snp_list = [ 
            {"chrom": "1", "start": 10000000, "end": 10000001, "variant_genotype": "T/T"}, 
            {"chrom": "1", "start": 10000000, "end": 10000001, "variant_genotype": "C/T"} 
        ]
        for snp_dir in snp_list:
            self.retrieve_match_and_assert(vcf_reader=get_vcf_reader, snp_dir=snp_dir, exp_length=1, exp_bool=False)

    def test_multi_snp(self, get_vcf_reader):
        snp_dir = {"chrom": "1", "start": 10000004, "end": 10000005, "variant_genotype": "C/C"}
        self.retrieve_match_and_assert(vcf_reader=get_vcf_reader, snp_dir=snp_dir, exp_length=2, exp_bool=True)



class TestIndelGenotypes():
    def retrieve_match_and_assert(self, vcf_reader, snp_dir, exp_length, exp_bool):
        matches = retrieve_match_all_records(vcf_reader=vcf_reader, snp=snp_dir)
        assert(len(matches) == exp_length)
        assert(all(matches) == exp_bool)

    def test_single_del_ref_ref_match(self, get_vcf_reader):
        variant_dir = {"chrom": "2", "start": 20000000, "end": 20000001, "variant_genotype": "TA/TA"}
        self.retrieve_match_and_assert(vcf_reader=get_vcf_reader, snp_dir=variant_dir, exp_length=1, exp_bool=True)

    def test_single_del_ref_alt_match(self, get_vcf_reader):
        variant_dir = {"chrom": "2", "start": 20000002, "end": 20000003, "variant_genotype": "TA/T"}
        self.retrieve_match_and_assert(vcf_reader=get_vcf_reader, snp_dir=variant_dir, exp_length=1, exp_bool=True)

    def test_single_del_alt_ref_match(self, get_vcf_reader):
        variant_dir = {"chrom": "2", "start": 20000004, "end": 20000005, "variant_genotype": "TA/T"}
        self.retrieve_match_and_assert(vcf_reader=get_vcf_reader, snp_dir=variant_dir, exp_length=1, exp_bool=True)

    def test_single_del_alt_alt_match(self, get_vcf_reader):
        variant_dir = {"chrom": "2", "start": 20000006, "end": 20000007, "variant_genotype": "T/T"}
        self.retrieve_match_and_assert(vcf_reader=get_vcf_reader, snp_dir=variant_dir, exp_length=1, exp_bool=True)

    def test_single_del_mismatch(self, get_vcf_reader):
        variant_list = [
            {"chrom": "2", "start": 20000000, "end": 20000001, "variant_genotype": "TA/T"},
            {"chrom": "2", "start": 20000000, "end": 20000001, "variant_genotype": "T/T"},
        ]
        for variant_dir in variant_list:
            self.retrieve_match_and_assert(vcf_reader=get_vcf_reader, snp_dir=variant_dir, exp_length=1, exp_bool=False)

    def test_single_ins_ref_ref_match(self, get_vcf_reader):
        variant_dir = {"chrom": "3", "start": 30000000, "end": 30000001, "variant_genotype": "G/G"}
        self.retrieve_match_and_assert(vcf_reader=get_vcf_reader, snp_dir=variant_dir, exp_length=1, exp_bool=True)

    def test_single_ins_ref_alt_match(self, get_vcf_reader):
        variant_dir = {"chrom": "3", "start": 30000002, "end": 30000003, "variant_genotype": "G/GA"}
        self.retrieve_match_and_assert(vcf_reader=get_vcf_reader, snp_dir=variant_dir, exp_length=1, exp_bool=True)

    def test_single_ins_alt_ref_match(self, get_vcf_reader):
        variant_dir = {"chrom": "3", "start": 30000004, "end": 30000005, "variant_genotype": "G/GA"}
        self.retrieve_match_and_assert(vcf_reader=get_vcf_reader, snp_dir=variant_dir, exp_length=1, exp_bool=True)

    def test_single_ins_alt_alt_match(self, get_vcf_reader):
        variant_dir = {"chrom": "3", "start": 30000006, "end": 30000007, "variant_genotype": "GA/GA"}
        self.retrieve_match_and_assert(vcf_reader=get_vcf_reader, snp_dir=variant_dir, exp_length=1, exp_bool=True)

    def test_single_ins_mismatch(self, get_vcf_reader):
        variant_list = [
            {"chrom": "3", "start": 30000000, "end": 30000001, "variant_genotype": "G/GA"},
            {"chrom": "3", "start": 30000000, "end": 30000001, "variant_genotype": "GA/GA"},
        ]
        for variant_dir in variant_list:
            self.retrieve_match_and_assert(vcf_reader=get_vcf_reader, snp_dir=variant_dir, exp_length=1, exp_bool=False)

    def test_multi_ins(self, get_vcf_reader):
        variant_dir = {"chrom": "3", "start": 30000008, "end": 30000009, "variant_genotype": "C/C"}
        self.retrieve_match_and_assert(vcf_reader=get_vcf_reader, snp_dir=variant_dir, exp_length=2, exp_bool=True)
