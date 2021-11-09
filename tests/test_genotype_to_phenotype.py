import pytest
import vcf as pyvcf

import assets.variant_genotype_to_phenotype as gt_to_pt

@pytest.fixture(scope="module", autouse=True)
def get_vcf_reader(setup_and_get_test_path):
    vcf_reader = pyvcf.Reader(filename=setup_and_get_test_path + "/fake_000000000000_R00C00.vcf.gz")
    return(vcf_reader)


class TestInput():
    def test_parser_required_args(self, setup_and_get_test_path):
        parser = gt_to_pt.parse_arguments_and_check(
            args_in=[setup_and_get_test_path + "fake_000000000000_R00C00.vcf.gz", "fake_000000000000_R00C00", "./references/sites_of_interest.yaml"])
        assert parser

    def test_parser_required_args_vcf(self, setup_and_get_test_path):
        parser = gt_to_pt.parse_arguments_and_check(
            args_in=[setup_and_get_test_path + "fake_000000000001_R01C01.vcf", "fake_000000000000_R00C00", "./references/sites_of_interest.yaml"])
        assert parser
    
    def test_parser_non_existing_input(self, setup_and_get_test_path):
        with pytest.raises(FileNotFoundError):
            parser = gt_to_pt.parse_arguments_and_check(
                args_in=[setup_and_get_test_path + "non_existing.vcf.gz", "fake_000000000000_R00C00", "./references/sites_of_interest.yaml"])

    def test_parser_non_existing_yaml(self, setup_and_get_test_path):
        with pytest.raises(FileNotFoundError):
            parser = gt_to_pt.parse_arguments_and_check(
                args_in=[setup_and_get_test_path + "fake_000000000000_R00C00.vcf.gz", "fake_000000000000_R00C00", "./references/non_existing.yaml"])

    def test_parser_unsupported_table(self, setup_and_get_test_path):
        with pytest.raises(TypeError):
            parser = gt_to_pt.parse_arguments_and_check(
                args_in=[setup_and_get_test_path + "fake_000000000000_R00C00.vcf.gz", "fake_000000000000_R00C00", "./references/sites_of_interest.json"])

    def test_parser_empty_input(self, setup_and_get_test_path):
        with pytest.raises(ValueError) as no_records_error:
            gt_to_pt.read_vcf(setup_and_get_test_path + "empty.vcf.gz")
        assert "empty" in str(no_records_error.value)

    def test_parser_empty_table(self, setup_and_get_test_path):
        with pytest.raises(ValueError) as empty_error:
            gt_to_pt.read_table(setup_and_get_test_path + "empty.yaml")
        assert "empty" in str(empty_error.value)

    def test_parser_no_records_input(self, setup_and_get_test_path):
        with pytest.raises(ValueError) as no_records_error:
            gt_to_pt.read_vcf(setup_and_get_test_path + "no_records.vcf.gz")
        assert "no records" in str(no_records_error.value)

    def test_parser_output_path(self, setup_and_get_test_path):
        parser = gt_to_pt.parse_arguments_and_check(
            args_in=[setup_and_get_test_path + "fake_000000000000_R00C00.vcf.gz", "fake_000000000000_R00C00", "./references/sites_of_interest.yaml", "--output_path", "./test_output/"])
        assert parser
    
    def test_parser_non_existing_output_path(self, setup_and_get_test_path):
        with pytest.raises(FileNotFoundError):
            parser = gt_to_pt.parse_arguments_and_check(
                args_in=[setup_and_get_test_path + "fake_000000000000_R00C00.vcf.gz", "fake_000000000000_R00C00", "./references/sites_of_interest.yaml", "--output_path", "./fake_dir/"])

    def test_parser_output_prefix(self, setup_and_get_test_path):
        parser = gt_to_pt.parse_arguments_and_check(
            args_in=[setup_and_get_test_path + "fake_000000000000_R00C00.vcf.gz", "fake_000000000000_R00C00", "./references/sites_of_interest.yaml", "--output_prefix", "test_prefix"])
        assert parser    


class TestSnpGenotypes():
    def retrieve_match_and_assert(self, vcf_reader, snp_dir, exp_length, exp_bool):
        matches = gt_to_pt.retrieve_match_all_records(vcf_reader=vcf_reader, snp=snp_dir)
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
        matches = gt_to_pt.retrieve_match_all_records(vcf_reader=vcf_reader, snp=snp_dir)
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
