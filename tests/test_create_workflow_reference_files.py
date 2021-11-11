#! venv/bin/python
# standard libraries alphabetic order of main package.
import configparser
import shutil

# third party libraries alphabetic order of main package.
import pytest

# local libraries alphabetic order of main package.
import assets.create_workflow_reference_files as create_ref

@pytest.fixture(scope="module", autouse=True)
def setup_and_get_test_path(tmp_path_factory):
    test_tmp_path = str(tmp_path_factory.mktemp("data")) + "/"

    # create empty config
    open(test_tmp_path + "empty.ini", "a").close()
    open(test_tmp_path + "empty_translation.csv", "a").close()
    # create config missing key
    config = configparser.ConfigParser()
    config['DEFAULT'] = {
        'species': 'homo sapiens',
        'translation_table': 'fake.csv',
    }
    with open(test_tmp_path + 'missing_key.ini', 'w') as config_missing_key:
        config.write(config_missing_key)

    # create config non existing translation file
    config['DEFAULT']['ensembl_url'] = "url"
    with open(test_tmp_path + 'non_existing_translation_file.ini', 'w') as config_no_tf:
        config.write(config_no_tf)
    
    # create config empty translation file
    config['DEFAULT']['translation_table'] = test_tmp_path + "empty_translation.csv"
    with open(test_tmp_path + 'empty_translation_file.ini', 'w') as config_empty_tf:
        config.write(config_empty_tf)

    shutil.copy("./assets/create_workflow_reference_files.ini", test_tmp_path)
    return test_tmp_path

# TODO implement these tests.
'''
rsID not exist in ensembl
rsID not linked to gene in ensembl
rsID and/or location linked to multiple genes -> not unique for single genotype gene.
variant ensembl query no result
variant with multiple alternative alleles
'''

class TestCreateRefsParser():
    def test_parser_optional_config_section(self):
        parser = create_ref.parse_arguments_and_check(
            args_in=["--config_section", "fake"]
        )
        assert parser
    
    def test_parser_optional_bed(self):
        parser = create_ref.parse_arguments_and_check(
            args_in=["--bed", "./references/sites_of_interest.bed"]
        )
        assert parser

    def test_parser_optional_output_path(self):
        parser = create_ref.parse_arguments_and_check(
            args_in=["--output_path", "./references/"]
        )
        assert parser

    def test_parser_optional_output_prefix(self):
        parser = create_ref.parse_arguments_and_check(
            args_in=["--output_prefix", "fake_prefix"]
        )
        assert parser

    def test_parser_optional_yaml(self):
        parser = create_ref.parse_arguments_and_check(
            args_in=["--yaml", "./references/sites_of_interest.yaml"]
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
                args_in=["--output_path", "./non_existing_path/"]
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

    def test_config_section_not_exists(self, setup_and_get_test_path):
        with pytest.raises(KeyError):
            config_section = create_ref.read_config_section_and_check(
                section="non_existing_section", config_file=setup_and_get_test_path + "create_workflow_reference_files.ini"
            )

    def test_config_key_not_exists(self, setup_and_get_test_path):
        with pytest.raises(KeyError):
            config_section = create_ref.read_config_section_and_check(
                section="DEFAULT", config_file=setup_and_get_test_path + "missing_key.ini"
            )

    def test_config_file_empty(self, setup_and_get_test_path):
        with pytest.raises(OSError) as empty_error:
            config_section = create_ref.read_config_section_and_check(
                section="DEFAULT", config_file=setup_and_get_test_path + "empty.ini"
            )
        assert "empty" in str(empty_error.value)
    
    def test_config_translation_empty(self, setup_and_get_test_path):
        with pytest.raises(OSError) as empty_error:
            config_section = create_ref.read_config_section_and_check(
                section="DEFAULT", config_file=setup_and_get_test_path + "empty_translation_file.ini"
            )
        assert "empty" in str(empty_error.value)

    def test_config_translation_not_exists(self, setup_and_get_test_path):
        with pytest.raises(FileNotFoundError):
            config_section = create_ref.read_config_section_and_check(
                section="DEFAULT", config_file=setup_and_get_test_path + "non_existing_translation_file.ini"
            )


# TODO implement these tests in TestCreateRefsTranslationFile.
'''
not all rename columns present.
genotype col does not contain : or /
genotype id not present in all rows.
gene_and_rs_id values not contain _
'''
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

