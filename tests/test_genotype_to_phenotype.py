#! venv/bin/python
# standard libraries alphabetic order of main package.
from pathlib import Path, PurePath
import shutil

# third party libraries alphabetic order of main package.
import pandas as pd
import pysam
import pytest
import vcf as pyvcf

# local libraries alphabetic order of main package.
import assets.variant_genotype_to_phenotype as gt_to_pt

@pytest.fixture(scope="module", autouse=True)
def setup_and_get_test_path(tmp_path_factory):
    test_tmp_path = str(tmp_path_factory.mktemp("tests")) + "/"
    shutil.rmtree(test_tmp_path)
    shutil.copytree("tests/test_data/", test_tmp_path)
    test_tmp_path_vcf = test_tmp_path + "/vcf_files/"
    # create empty files
    open(str(test_tmp_path) + "/empty.yaml", "a").close()
    open(str(test_tmp_path) + "/empty.json", "a").close()
    for vcf_file in Path(test_tmp_path + "vcf_files").glob("**/*.vcf"):
        basename = PurePath(vcf_file).name
        tmp_vcf_file = str(test_tmp_path_vcf) + "/" + basename
        pysam.tabix_compress(tmp_vcf_file, tmp_vcf_file + ".gz")
        pysam.tabix_index(tmp_vcf_file + ".gz", preset="vcf")
    shutil.copy(str(test_tmp_path_vcf) + "/" + "sample.vcf",
                str(test_tmp_path_vcf) + "/" + "sample_copy.vcf")
    return str(test_tmp_path) + "/"


@pytest.fixture(scope="module", autouse=True)
def get_vcf_reader(setup_and_get_test_path):
    vcf_reader = pyvcf.Reader(filename=setup_and_get_test_path + "/vcf_files/sample.vcf.gz")
    return(vcf_reader)


@pytest.fixture(scope="module")
def get_tmp_output_path(tmp_path_factory):
    output_path = str(tmp_path_factory.mktemp("output")) + "/"
    # open(str(output_path) + "existing_output.bed", "a").close()
    return output_path


# TODO: implement tests:
## pos 1 del, pos 2 ins, only match on expected genotype measurement type. (fetch will pick up both.)

class TestGtToPtInputs():
    def test_parser_required_args(self, setup_and_get_test_path):
        parser = gt_to_pt.parse_arguments_and_check(
            args_in=[setup_and_get_test_path + "/vcf_files/sample.vcf.gz", "sample", "./references/sites_of_interest.yaml"])
        assert parser

    def test_parser_required_args_vcf(self, setup_and_get_test_path):
        parser = gt_to_pt.parse_arguments_and_check(
            args_in=[setup_and_get_test_path + "/vcf_files/sample_copy.vcf", "sample", "./references/sites_of_interest.yaml"])
        assert parser
    
    def test_parser_non_existing_input(self, setup_and_get_test_path):
        with pytest.raises(FileNotFoundError):
            parser = gt_to_pt.parse_arguments_and_check(
                args_in=[setup_and_get_test_path + "non_existing.vcf.gz", "sample", "./references/sites_of_interest.yaml"])

    def test_parser_non_existing_yaml(self, setup_and_get_test_path):
        with pytest.raises(FileNotFoundError):
            parser = gt_to_pt.parse_arguments_and_check(
                args_in=[setup_and_get_test_path + "/vcf_files/sample.vcf.gz", "sample", "./references/non_existing.yaml"])

    def test_parser_unsupported_table(self, setup_and_get_test_path):
        with pytest.raises(TypeError):
            parser = gt_to_pt.parse_arguments_and_check(
                args_in=[setup_and_get_test_path + "/vcf_files/sample.vcf.gz", "sample", "./references/sites_of_interest.json"])

    def test_parser_empty_input(self, setup_and_get_test_path):
        with pytest.raises(ValueError) as no_records_error:
            gt_to_pt.read_vcf(setup_and_get_test_path + "/vcf_files/empty.vcf.gz")
        assert "empty" in str(no_records_error.value)

    def test_parser_empty_table(self, setup_and_get_test_path):
        with pytest.raises(ValueError) as empty_error:
            gt_to_pt.read_yaml(setup_and_get_test_path + "empty.yaml")
        assert "empty" in str(empty_error.value)

    def test_parser_no_records_input(self, setup_and_get_test_path):
        with pytest.raises(ValueError) as no_records_error:
            gt_to_pt.read_vcf(setup_and_get_test_path + "/vcf_files/no_records.vcf.gz")
        assert "no records" in str(no_records_error.value)

    def test_parser_output_path(self, setup_and_get_test_path):
        parser = gt_to_pt.parse_arguments_and_check(
            args_in=[setup_and_get_test_path + "/vcf_files/sample.vcf.gz", "sample", "./references/sites_of_interest.yaml", "--output_path", "./test_output/"])
        assert parser
    
    def test_parser_non_existing_output_path(self, setup_and_get_test_path):
        with pytest.raises(FileNotFoundError):
            parser = gt_to_pt.parse_arguments_and_check(
                args_in=[setup_and_get_test_path + "/vcf_files/sample.vcf.gz", "sample", "./references/sites_of_interest.yaml", "--output_path", "./fake_dir/"])

    def test_parser_output_prefix(self, setup_and_get_test_path):
        parser = gt_to_pt.parse_arguments_and_check(
            args_in=[setup_and_get_test_path + "/vcf_files/sample.vcf.gz", "sample", "./references/sites_of_interest.yaml", "--output_prefix", "test_prefix"])
        assert parser    

    def test_check_reqs_keys(self):
        gt_to_pt.check_required_keys(snp={"chrom": 0, "start": 0, "end": 1, "variant_genotype": "A/A", "name": "fake_id"})
        assert True

    def test_check_reqs_keys_missing_keys(self):
        dir_snp_list = [
            {"start": 0, "end": 1, "variant_genotype": "A/A", "name": "fake_id"}, # missing chrom
            {"chrom": 0, "end": 1, "variant_genotype": "A/A", "name": "fake_id"},  # missing start
            {"chrom": 0, "start": 0, "variant_genotype": "A/A", "name": "fake_id"}, # missing end
            {"chrom": 0, "start": 0, "end": 1, "name": "fake_id"}, # missing variant_genotype
            {"chrom": 0, "start": 0, "end": 1, "variant_genotype": "A/A"}, # missing variant_genotype
        ]
        for snp_dir in dir_snp_list:
            with pytest.raises(ValueError) as missing_key_error:
                gt_to_pt.check_required_keys(snp=snp_dir)
            assert "missing" in str(missing_key_error.value).lower()


class TestGtToPtSnpGenotypes():
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


class TestGtToPtIndelGenotypes():
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


class TestGtToPtMatchSnpGenotype():
    def test_match_snp_genotype_no_record(self, setup_and_get_test_path):
        vcf_reader = pyvcf.Reader(filename=setup_and_get_test_path + "/vcf_files/sample.vcf.gz")
        genotypes = gt_to_pt.read_yaml(
            translation_file=setup_and_get_test_path + "translation_yaml/gt_missing_rs.yaml"
        ).values()
        with (
            pytest.warns(UserWarning, match="fetch has no records"), 
            pytest.warns(UserWarning, match="No genotype match found."),
        ):
            matched_gt = gt_to_pt.retrieve_match_snp_genotype(vcf_reader=vcf_reader, genotypes=genotypes)
        assert not matched_gt

    def test_match_snp_genotype_gt_single_rs(self, setup_and_get_test_path):
        vcf_reader = pyvcf.Reader(filename=setup_and_get_test_path + "/vcf_files/sample.vcf.gz")
        genotypes = gt_to_pt.read_yaml(
            translation_file=setup_and_get_test_path + "translation_yaml/gt_single_rs.yaml"
        ).values()
        matched_gt = gt_to_pt.retrieve_match_snp_genotype(vcf_reader=vcf_reader, genotypes=genotypes)
        assert '0' in matched_gt

    def test_match_snp_genotype_gt_multi_rs(self, setup_and_get_test_path):
        vcf_reader = pyvcf.Reader(filename=setup_and_get_test_path + "/vcf_files/sample.vcf.gz")
        genotypes = gt_to_pt.read_yaml(
            translation_file=setup_and_get_test_path + "translation_yaml/gt_multi_rs.yaml"
        ).values()
        matched_gt = gt_to_pt.retrieve_match_snp_genotype(vcf_reader=vcf_reader, genotypes=genotypes)
        assert '0' in matched_gt
        assert matched_gt.get('0') == [True, True]
        assert '1' in matched_gt
        assert matched_gt.get('1') == [True, False]


class TestGtToPtWriteMatchedPhenotype():
    def count_lines(self, filename):
        with open(filename) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def test_write_matched_pt_check_format_output(self, setup_and_get_test_path, get_tmp_output_path):
        matched_gt = {'0': [True]}
        translation = gt_to_pt.read_yaml(
            translation_file=setup_and_get_test_path + "/translation_yaml/gt_single_rs.yaml"
        )
        gt_to_pt.write_matched_phenotype(
            genotype_match_per_snp=matched_gt, 
            translation_table=translation, 
            output_path=get_tmp_output_path,
            output_prefix="test_output_sample", 
            sample="sample"
        )
        assert Path(get_tmp_output_path + "/test_output_sample.txt").is_file()
        data = pd.read_csv(get_tmp_output_path + "/test_output_sample.txt", sep="\t")
        expected_columns = ["sample", "genotype_match", "phenotype_id", "phenotype_match"]
        assert len(set(expected_columns).intersection(data)) == len(expected_columns)
        assert len(data.index) == 1 # test single match aka single row.
        assert data.loc[0,'genotype_match'] == "fakegene:wildtype/wildtype"
        assert data.loc[0, "phenotype_id"] == 0
        assert data.loc[0, "phenotype_match"] == "POOR METABOLIZER"

    def test_write_matched_pt_single_match(self, setup_and_get_test_path, get_tmp_output_path):
        matched_gt = {'0': [True, True]}
        translation = gt_to_pt.read_yaml(
            translation_file=setup_and_get_test_path + "/translation_yaml/gt_multi_rs.yaml"
        )
        gt_to_pt.write_matched_phenotype(
            genotype_match_per_snp=matched_gt,
            translation_table=translation,
            output_path=get_tmp_output_path,
            output_prefix="test_output_sample",
            sample="sample"
        )
        n_lines = self.count_lines(filename=get_tmp_output_path + "/test_output_sample.txt")
        assert n_lines == 2

    def test_write_matched_pt_two_gt_single_match(self, setup_and_get_test_path, get_tmp_output_path):
        matched_gt = {'0': [True, True], '1': [True, False]}
        translation = gt_to_pt.read_yaml(
            translation_file=setup_and_get_test_path + "/translation_yaml/gt_multi_rs.yaml"
        )
        gt_to_pt.write_matched_phenotype(
            genotype_match_per_snp=matched_gt,
            translation_table=translation,
            output_path=get_tmp_output_path,
            output_prefix="test_output_sample",
            sample="sample"
        )
        n_lines = self.count_lines(filename=get_tmp_output_path + "/test_output_sample.txt")
        assert n_lines == 2

    def test_write_matched_pt_multi_match(self, setup_and_get_test_path, get_tmp_output_path):
        matched_gt = {'0': [True, True], '1': [True, True]}
        translation = gt_to_pt.read_yaml(
            translation_file=setup_and_get_test_path + "/translation_yaml/gt_multi_rs.yaml"
        )
        gt_to_pt.write_matched_phenotype(
            genotype_match_per_snp=matched_gt,
            translation_table=translation,
            output_path=get_tmp_output_path,
            output_prefix="test_output_sample",
            sample="sample"
        )
        n_lines = self.count_lines(filename=get_tmp_output_path + "/test_output_sample.txt")
        assert n_lines == 3

