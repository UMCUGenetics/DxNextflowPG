#! venv/bin/python
# standard libraries alphabetic order of main package.
import configparser
import csv
import shutil
import pandas as pd

# third party libraries alphabetic order of main package.
import requests
import pytest

# local libraries alphabetic order of main package.
import assets.create_workflow_reference_files as create_ref


# TODO implement these tests.
"""
rsID and/or location linked to multiple genes -> not unique for single genotype gene.
variant with multiple alternative alleles ensembl
variant with multiple alternative alleles translation table
"""


@pytest.fixture(scope="module")
def create_test_files(tmp_path_factory):
    test_path = str(tmp_path_factory.mktemp("tests")) + "/"
    shutil.rmtree(test_path)
    shutil.copytree("tests/test_data/", test_path)
    return test_path

class TestCreateRefsParser():
    def test_parser_optional_config_section(self):
        parser = create_ref.parse_arguments_and_check(
            args_in=["--config_section", "fake"]
        )
        assert parser
    
    def test_parser_optional_bed(self):
        parser = create_ref.parse_arguments_and_check(
            args_in=["--bed", "references/sites_of_interest.bed"]
        )
        assert parser

    def test_parser_optional_output_path(self):
        parser = create_ref.parse_arguments_and_check(
            args_in=["--output_path", "references/"]
        )
        assert parser

    def test_parser_optional_output_prefix(self):
        parser = create_ref.parse_arguments_and_check(
            args_in=["--output_prefix", "fake_prefix"]
        )
        assert parser

    def test_parser_optional_yaml(self):
        parser = create_ref.parse_arguments_and_check(
            args_in=["--yaml", "references/sites_of_interest.yaml"]
        )
        assert parser

    def test_parser_optional_bed_not_exist(self):
        with pytest.raises(FileNotFoundError):
            parser = create_ref.parse_arguments_and_check(
                args_in=["--bed", "non_existing_file.bed"]
            )

    def test_parser_optional_output_path_not_exist(self):
        with pytest.raises(FileNotFoundError):
            parser = create_ref.parse_arguments_and_check(
                args_in=["--output_path", "non_existing_path/"]
            )

    def test_parser_optional_yaml_not_exist(self):
        with pytest.raises(FileNotFoundError):
            parser = create_ref.parse_arguments_and_check(
                args_in=["--yaml", "non_existing_file.yaml"]
            )


class TestCreateRefsConfig():
    def test_config_file_not_exists(self):
        with pytest.raises(FileNotFoundError):
            config_section = create_ref.read_config_section_and_check(
                section="DEFAULT", config_file="non_existing_file.ini"
            )

    def test_config_section_not_exists(self, create_test_files):
        with pytest.raises(KeyError):
            config_section = create_ref.read_config_section_and_check(
                section="non_existing_section", config_file=create_test_files + "/ini_config/correct_config.ini"
            )

    def test_config_key_not_exists(self, create_test_files):
        with pytest.raises(KeyError):
            config_section = create_ref.read_config_section_and_check(
                section="DEFAULT", config_file=create_test_files + "/ini_config/missing_key.ini"
            )

    def test_config_file_empty(self, create_test_files):
        with pytest.raises(OSError) as empty_error:
            config_section = create_ref.read_config_section_and_check(
                section="DEFAULT", config_file=create_test_files + "/ini_config/empty.ini"
            )
        assert "empty" in str(empty_error.value)
    
    def test_config_translation_empty(self, create_test_files):
        with pytest.raises(OSError) as empty_error:
            config_section = create_ref.read_config_section_and_check(
                section="DEFAULT", config_file=create_test_files + "/ini_config/empty_translation_file.ini"
            )
        assert "empty" in str(empty_error.value)

    def test_config_translation_not_exists(self, create_test_files):
        with pytest.raises(FileNotFoundError):
            config_section = create_ref.read_config_section_and_check(
                section="DEFAULT", config_file=create_test_files + "/ini_config/non_existing_translation_file.ini"
            )


class TestCreateRefsTranslationFile():
    pass


# TODO implement these tests in TestCreateRefsEnsembl.
'''
connection ensembl cannot be established (typo url)
'''
class TestCreateRefsEnsembl():
    pass


# TODO implement these tests in TestCreateRefsInvalidGenes
'''
genotype gene with multiple rs IDs, but from different chromosomes.
genotype gene not match with ensembl retrieved gene
'''
class TestCreateRefsInvalidGenes():
    pass


# TODO implement these tests in TestCreateRefsFilterGenotypes
'''
if no dict_rename_tf_cols, succeed.
if no filter_rs_id, succeed.
if no filter_gene_regex, succeed
filter of rs ID should remove all related genotypes.
filter of gene should remove all related genotypes.
no genotypes left after filter rs_ids and genes --> raise ERROR
'''
class TestCreateRefsFilterGenotypes():
    pass


# TODO implement these tests in TestCreateRefsVariantGenotypes
'''
Indel notation not using ., but another special character.
Indel notation not using ., already 'correct'.

variant genotypes 1 bp on reverse strand are changed to forward notation.
- REF/REF
- REF/ALT
- ALT/REF
- ALT/ALT
variant genotypes > 1bp on reverse strand are changed to forward notation.
- REF/REF
- REF/ALT
- ALT/REF
- ALT/ALT
variants genotypes indels notation standardized.
- REF/REF
- REF/ALT
- ALT/REF
- ALT/ALT
'''
class TestCreateRefsVariantGenotypes():
    pass


# TODO implement these tests in TestCreateRefsCompareFiles
'''
compare files, not supported data type.
compare files, dict yaml, only order different dict
compare files, dict yaml, different variant 
compare files, dict yaml, removed genotype
compare files, dict yaml, removed random field.
compare files, bed, removed rs site.
compare files, bed, synonym rs id and thus name.
'''
class TestCreateRefsCompareFiles():
    pass


# TODO implement these tests in TestCreateRefsYaml
'''
yaml write, file exists already. - -> ERROR, param 'force'?
yaml write, file exists already, use force - -> succeed.
'''
class TestCreateRefsYaml():
    pass


# TODO implement these tests in TestCreateRefsBed
'''
write bedfile to non existing file - -> succeed
write bedfile, file exists already - -> ERROR, param 'force'?
'''
class TestCreateRefsBed():
    pass


# TODO implement these tests in TestCreateRefsOutput
'''
output file starts with output_prefix.
output file saved in output path.
'''
class TestCreateRefsOutput():
    pass

